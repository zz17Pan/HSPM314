classdef MultiModalSensing < handle
    properties
        cfg                 % ϵͳ���ò���
        dict               % �����ֵ�
        sensing_params     % ��֪����
        sensor_state      % ��ǰ��֪״̬
        valid_frames = 0   % ��Ч֡����
        snr_history = []   % SNR��ʷ��¼
        error_history = [] % �����ʷ��¼
    end
    
    methods
        % ���캯��
        function obj = MultiModalSensing()
            try
                obj.cfg = Config;
                
                % ���ø�֪����
                obj.sensing_params = struct(...
                    'fmcw_bandwidth', 1e9, ...    % 1GHz����
                    'sweep_time', 1e-3, ...       % 1msɨ��ʱ��
                    'N_samples', 1024, ...        % ��������
                    'N_range_fft', 2048, ...      % ����FFT����
                    'N_doppler_fft', 256, ...     % ������FFT����
                    'dict_size', [180, 180, 50],... % �ֵ�ά��
                    'omp_sparsity', 3, ...        % OMPϡ���
                    'min_snr_db', 5);             % ��СSNR��ֵ
                
                fprintf('ϵͳ����: Nx=%d, Nz=%d, ��������=%d\n', ...
                    obj.cfg.Nx, obj.cfg.Nz, obj.cfg.Nx * obj.cfg.Nz);
                
                fprintf('��ʼ���������ֵ�...\n');
                obj.dict = obj.generate_dictionary();
                fprintf('�ֵ��������\n');
                
            catch ME
                fprintf('��ʼ��ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % ��ģ̬��֪������
        function result = perform_multimodal_sensing(obj, tx_pos, rx_pos, t)
            try
                fprintf('\n����� %d ֡ (t = %.3f s)\n', obj.valid_frames + 1, t);
                
                % 1. �Ƕ����֪
                angle_result = obj.perform_angle_sensing(tx_pos, rx_pos, t);
                
                % ��֤���������
                required_fields = {'position', 'theta', 'phi', 'range', 'received_signal'};
                for i = 1:length(required_fields)
                    if ~isfield(angle_result, required_fields{i})
                        error('�Ƕȸ�֪���ȱ�ٱ�Ҫ�ֶ�: %s', required_fields{i});
                    end
                end
                
                % 2. ����-�����ո�֪
                range_doppler = obj.perform_range_doppler_sensing(tx_pos, rx_pos);
                
                % 3. �������֪
                polarization = obj.perform_polarization_sensing(tx_pos, rx_pos);
                
                % 4. �����ں�
                fusion_result = obj.feature_fusion(angle_result, range_doppler, polarization);
                
                % �����������
                result = struct(...
                    'timestamp', t, ...
                    'position', fusion_result.position, ...
                    'velocity', range_doppler.velocity, ...
                    'theta', angle_result.theta, ...
                    'phi', angle_result.phi, ...
                    'range', angle_result.range, ...
                    'polarization', polarization, ...
                    'confidence', fusion_result.confidence, ...
                    'received_signal', angle_result.received_signal, ...
                    'range_profile', range_doppler.range_profile, ...
                    'doppler_profile', range_doppler.doppler_profile);
                
                obj.valid_frames = obj.valid_frames + 1;
                obj.sensor_state = result;
                obj.update_performance_metrics(result);
                
            catch ME
                fprintf('��ģ̬��֪����ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % �Ƕ����֪
        function result = perform_angle_sensing(obj, tx_pos, rx_pos, t)
            try
                % 1. ��ȡ��������
                [~, sensing_subarrays] = ResourceAllocation.allocate_subarrays();
                
                % 2. ���ɵ�Ƶ�ź�
                N_tx = size(tx_pos, 1);
                pilot = obj.generate_pilot_signal(N_tx);
                
                % 3. �����ŵ�����
                H = ChannelHSPM.generate_channel(tx_pos, rx_pos, t);
                
                % 4. �����źŴ���
                Y = H * pilot;
                
                % 5. ά��ƥ�䴦��
                N_ant = obj.cfg.Nx * obj.cfg.Nz;
                if size(Y,1) > N_ant
                    Y = Y(1:N_ant, :);
                end
                
                % 6. SNR����
                snr = obj.estimate_snr(Y);
                fprintf('����SNR: %.2f dB\n', snr);
                
                % 7. ��������
                if snr >= obj.sensing_params.min_snr_db
                    [theta, phi, r] = obj.estimate_parameters(Y, obj.dict);
                else
                    fprintf('SNR������ֵ��ʹ����һ֡���\n');
                    last_est = obj.get_last_estimate();
                    theta = last_est.theta;
                    phi = last_est.phi;
                    r = last_est.range;
                end
                
                % 8. ������
                result = struct(...
                    'position', [r*cos(phi)*cos(theta), r*cos(phi)*sin(theta), r*sin(phi)], ...
                    'theta', theta, ...
                    'phi', phi, ...
                    'range', r, ...
                    'velocity', 0, ...
                    'received_signal', Y);
                
                fprintf('�Ƕȹ��ƽ����\n');
                fprintf('��λ��: %.2f��\n', rad2deg(theta));
                fprintf('������: %.2f��\n', rad2deg(phi));
                fprintf('����: %.2f m\n', r);
                
            catch ME
                fprintf('�Ƕȸ�֪ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % ��������ո�֪
        function result = perform_range_doppler_sensing(obj, tx_pos, rx_pos)
            try
                % 1. FMCW��������
                params = obj.sensing_params;
                N_samples = params.N_samples;
                Ts = params.sweep_time / N_samples;
                t_samples = (0:N_samples-1) * Ts;
                
                % 2. ����FMCW�ź�
                fmcw = obj.generate_fmcw_signal(t_samples, params.fmcw_bandwidth, params.sweep_time);
                
                % 3. �����ŵ�����
                H = ChannelHSPM.generate_channel(tx_pos, rx_pos, 0);
                
                % 4. �����źŴ���
                tx_signal = repmat(fmcw(:), 1, size(H, 2))';
                rx_signal = H * tx_signal;
                
                % 5. ��ʼ������Ͷ���������
                N_rx = size(rx_signal, 1);
                range_profiles = zeros(N_rx, params.N_range_fft/2);
                doppler_profiles = zeros(N_rx, params.N_doppler_fft/2);
                
                % 6. ����ÿ�������ź�
                for n = 1:N_rx
                    [range_prof, doppler_prof] = obj.process_fmcw(...
                        rx_signal(n,:), params.N_range_fft, params.N_doppler_fft);
                    range_profiles(n,:) = range_prof;
                    doppler_profiles(n,:) = doppler_prof;
                end
                
                % 7. ������
                last_est = obj.get_last_estimate();
                result = struct(...
                    'range', last_est.range, ...
                    'velocity', 0, ...
                    'position', last_est.position, ...
                    'theta', last_est.theta, ...
                    'phi', last_est.phi, ...
                    'range_profile', mean(range_profiles, 1), ...
                    'doppler_profile', mean(doppler_profiles, 1), ...
                    'received_signal', rx_signal);
                
            catch ME
                fprintf('��������ո�֪ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % ����FMCW�ź�
        function signal = generate_fmcw_signal(~, t, bandwidth, sweep_time)
            % �����Ƶб��
            mu = bandwidth / sweep_time;
            
            % ���ɻ���FMCW�ź�
            signal = exp(1j*pi*mu*t.^2);
            
            % ��Ӵ������͹�һ��
            window = hamming(length(t));
            signal = signal(:) .* window;
            signal = signal / norm(signal);
        end

        % �������֪
        function result = perform_polarization_sensing(obj, tx_pos, rx_pos)
            try
                % ��ȡ�����ŵ�
                H = ChannelHSPM.generate_channel(tx_pos, rx_pos, 0);
                
                % ��ȡ��������
                N_ant = obj.cfg.Nx * obj.cfg.Nz;
                H_HH = H(1:N_ant/2, 1:N_ant/2);
                H_HV = H(1:N_ant/2, N_ant/2+1:end);
                H_VH = H(N_ant/2+1:end, 1:N_ant/2);
                H_VV = H(N_ant/2+1:end, N_ant/2+1:end);
                
                % ���켫������
                H_pol = [mean2(H_HH) mean2(H_HV); mean2(H_VH) mean2(H_VV)];
                
                % ���Ƽ�������
                [U, S, V] = svd(H_pol);
                pol_angle = angle(V(1,1));
                ellipticity = S(2,2) / S(1,1);
                
                result = struct(...
                    'pol_angle', pol_angle, ...
                    'ellipticity', ellipticity, ...
                    'received_signal', H);
                
            catch ME
                fprintf('������֪ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % �����ں�
        function result = feature_fusion(obj, angle_result, range_doppler, polarization)
            try
                % ������������Ŷ�
                conf_angle = obj.calculate_confidence(obj.estimate_snr(angle_result.received_signal));
                conf_range = obj.calculate_confidence(obj.estimate_snr(range_doppler.received_signal));
                conf_pol = obj.calculate_confidence(obj.estimate_snr(polarization.received_signal));
                
                % D-S֤�������ں�
                weights = [conf_angle, conf_range, conf_pol];
                fusion_coef = obj.ds_fusion(weights);
                
                % λ���ں�
                pos_angle = angle_result.position;
                pos_range = [range_doppler.range * cos(angle_result.theta), ...
                           range_doppler.range * sin(angle_result.theta), ...
                           range_doppler.range * sin(angle_result.phi)];
                           
                fused_position = fusion_coef * pos_angle + (1-fusion_coef) * pos_range;
                
                result = struct(...
                    'position', fused_position, ...
                    'confidence', fusion_coef);
                
            catch ME
                fprintf('�����ں�ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % �����ֵ�
        function dict = generate_dictionary(obj)
            try
                % 1. ��ȡϵͳ����
                Nx = obj.cfg.Nx;
                Nz = obj.cfg.Nz;
                N_ant = Nx * Nz;
                dict_size = obj.sensing_params.dict_size;
                total_atoms = dict_size(1) * dict_size(2) * dict_size(3);
                
                % 2. ��ʼ���ֵ����
                dict = zeros(N_ant, total_atoms);
                fprintf('��ʼ�����ֵ䣬ԭ������: %d\n', total_atoms);
                
                % 3. ���ɲ�������
                theta_grid = linspace(-pi/2, pi/2, dict_size(1));
                phi_grid = linspace(-pi/2, pi/2, dict_size(2));
                r_grid = linspace(1, 100, dict_size(3));
                
                % 4. �����ֵ�ԭ��
                idx = 1;
                progress_step = floor(total_atoms/10);
                
                for i_theta = 1:dict_size(1)
                    for i_phi = 1:dict_size(2)
                        for i_r = 1:dict_size(3)
                            dict(:,idx) = obj.generate_atom(...
                                theta_grid(i_theta), ...
                                phi_grid(i_phi), ...
                                r_grid(i_r));
                            
                            if mod(idx, progress_step) == 0
                                fprintf('�ֵ����ɽ���: %.1f%%\n', 100*idx/total_atoms);
                            end
                            idx = idx + 1;
                        end
                    end
                end
                
                % 5. �������͹�һ��
                fprintf('����������...\n');
                [dict, ~] = qr(dict, 0);
                for i = 1:size(dict, 2)
                    dict(:,i) = dict(:,i) / norm(dict(:,i));
                end
                
                fprintf('�ֵ�ά��: [%d, %d]\n', size(dict));
                
            catch ME
                fprintf('�ֵ�����ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
%         % �����ֵ�ԭ��
%         function atom = generate_atom(obj, theta, phi, r)
%             try
%                 % 1. ��ȡϵͳ����
%                 Nx = obj.cfg.Nx;
%                 Nz = obj.cfg.Nz;
%                 d = obj.cfg.d;
%                 lambda = obj.cfg.lambda;
%                 k = 2*pi/lambda;
%                                 % ���� generate_atom ����
%                 % 2. ��ʼ��ԭ������
%                 atom = zeros(Nx*Nz, 1);
%                 idx = 1;
%                 
%                 % 3. ���ɿռ���Ӧ
%                 for nz = 0:(Nz-1)
%                     for nx = 0:(Nx-1)
%                         % ������λ
%                         phase = k * d * (...
%                             nx * sin(theta) * cos(phi) + ...
%                             nz * sin(phi));
%                         
%                         % ��Ӿ���˥��
%                         amp = exp(-r/100);
%                         
%                         % ���ɸ���Ӧ
%                         atom(idx) = amp * exp(1j*phase);
%                         idx = idx + 1;
%                     end
%                 end
%                 
%                 % 4. ��һ��
%                 atom = atom / norm(atom);
%                 
%             catch ME
%                 fprintf('ԭ������ʧ��: %s\n', ME.message);
%                 rethrow(ME);
%             end
%         end
function atom = generate_atom(obj, theta, phi, r)
    % 1. ��ȡϵͳ����
    Nx = obj.cfg.Nx;
    Nz = obj.cfg.Nz;
    d = obj.cfg.d;
    lambda = obj.cfg.lambda;
    k = 2*pi/lambda;
    
    % 2. ���ɿռ�·�����ģ��
    path_loss = (lambda/(4*pi*r))^2;  % ����˥��
    amp = sqrt(path_loss);            % ����˥��
    
    % 3. ��ʼ��ԭ������
    atom = zeros(Nx*Nz, 1);
    idx = 1;
    
    % 4. ���ɿռ���Ӧ
    for nz = 0:(Nz-1)
        for nx = 0:(Nx-1)
            % ��λ���㣨������ά���Σ�
            phase = k * d * (...
                nx * sin(theta) * cos(phi) + ...
                nz * sin(phi)...
            );
            atom(idx) = amp * exp(1j*phase);
            idx = idx + 1;
        end
    end
    
    % 5. ��һ��
    atom = atom / norm(atom);
end
        % OMP�㷨ʵ��
        function x = omp(obj, dict, y, K)
            try
                % 1. ����Ԥ��������
                [M, N] = size(dict);
                y = y(:);  % ȷ��y��������
                
                if length(y) ~= M
                    error('ά�Ȳ�ƥ��: dict(%d��%d), y(%d��1)', M, N, length(y));
                end
                
                % 2. ��ʼ��
                x = zeros(N, 1);
                residual = y;
                support = [];
                D_s = [];
                x_s = [];
                
                % 3. �����ؽ�
                for k = 1:min(K, N)
                    % ���������
                    corr = abs(dict' * residual);
                    [~, idx] = max(corr);
                    
                    % ����֧�ּ�
                    if ~ismember(idx, support)
                        support = [support; idx];
                        
                        % ��С���˹���
                        D_s = dict(:, support);
                        if rank(D_s) < length(support)
                            support = support(1:end-1);
                            continue;
                        end
                        
                        x_s = pinv(D_s) * y;
                        
                        % ���²в�
                        residual = y - D_s * x_s;
                        
                        % �������
                        if norm(residual) < 1e-6
                            break;
                        end
                    end
                end
                
                % 4. �����
                x(support) = x_s;
                
            catch ME
                fprintf('OMP�㷨ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % ��������
        function [theta, phi, r] = estimate_parameters(obj, Y, dict)
            try
                % 1. �������ά��
                [M1, N1] = size(dict);
                [M2, N2] = size(Y);
                
                % 2. ȷ��Y��������
                if N2 > 1
                    Y = Y(:,1);
                end
                
                % 3. ��֤ά��ƥ��
                if M1 ~= M2
                    error('�ֵ�͹۲��ź�ά�Ȳ�ƥ�䣺dict(%d��%d), Y(%d��%d)', ...
                        M1, N1, M2, N2);
                end
                
                % 4. ִ��OMP�ؽ�
                x = obj.omp(dict, Y, obj.sensing_params.omp_sparsity);
                
                % 5. �ҵ���ǿ����
                [~, max_idx] = max(abs(x));
                
                % 6. �������ָ�����
                dict_size = obj.sensing_params.dict_size;
                [i_theta, i_phi, i_r] = ind2sub(dict_size, max_idx);
                
                % 7. ����ʵ�ʲ���ֵ
                theta_grid = linspace(-pi/2, pi/2, dict_size(1));
                phi_grid = linspace(-pi/2, pi/2, dict_size(2));
                r_grid = linspace(1, 100, dict_size(3));
                
                theta = theta_grid(i_theta);
                phi = phi_grid(i_phi);
                r = r_grid(i_r);
                
            catch ME
                fprintf('��������ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % ���ɵ�Ƶ�ź�
        function pilot = generate_pilot_signal(~, N_tx)
            try
                pilot = eye(N_tx);
                phases = exp(1j*2*pi*rand(N_tx, 1));
                pilot = pilot * diag(phases);
                pilot = pilot / sqrt(N_tx);
                
            catch ME
                fprintf('��Ƶ�ź�����ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % FMCW�źŴ���
        function [range_profile, doppler_profile] = process_fmcw(~, rx_signal, N_range_fft, N_doppler_fft)
            try
                % 1. ȷ���ź�Ϊ������
                rx_signal = rx_signal(:).';
                
                % 2. Ӧ�ô�����
                window = hamming(length(rx_signal));
                rx_signal_windowed = rx_signal .* window.';
                
                % 3. ����FFT����
                range_fft = fft(rx_signal_windowed, N_range_fft);
                range_profile = abs(range_fft(1:N_range_fft/2));
                
                % 4. ������FFT����
                doppler_fft = fft(rx_signal_windowed, N_doppler_fft);
                doppler_profile = abs(doppler_fft(1:N_doppler_fft/2));
                
                % 5. ��һ��
                range_profile = range_profile / max(abs(range_profile) + eps);
                doppler_profile = doppler_profile / max(abs(doppler_profile) + eps);
                
            catch ME
                fprintf('FMCW����ʧ��: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % �������Ŷ�
        function conf = calculate_confidence(~, SNR)
            conf = 1 - exp(-SNR/10);
        end
        
        % D-S֤�������ں�
        function result = ds_fusion(~, masses)
            result = prod(masses) / (prod(masses) + prod(1-masses));
        end
        
        % ��ȡ��һ֡���ƽ��
        function last_est = get_last_estimate(obj)
            try
                if obj.valid_frames > 0 && ~isempty(obj.sensor_state)
                    last_est = struct(...
                        'position', obj.sensor_state.position, ...
                        'theta', obj.sensor_state.theta, ...
                        'phi', obj.sensor_state.phi, ...
                        'range', obj.sensor_state.range, ...
                        'velocity', obj.sensor_state.velocity);
                else
                    last_est = struct(...
                        'position', [0, 0, 40], ...
                        'theta', 0, ...
                        'phi', pi/4, ...
                        'range', 40, ...
                        'velocity', 0);
                end
            catch ME
                fprintf('��ȡ��һ֡����ʧ�ܣ�ʹ��Ĭ��ֵ\n');
                last_est = struct(...
                    'position', [0, 0, 40], ...
                    'theta', 0, ...
                    'phi', pi/4, ...
                    'range', 40, ...
                    'velocity', 0);
            end
        end
        
        % ��������ָ��
        function update_performance_metrics(obj, result)
            try
                current_snr = obj.estimate_snr(result.received_signal);
                obj.snr_history(end+1) = current_snr;
                
                if isfield(result, 'true_position') && isfield(result, 'position')
                    error = norm(result.position - result.true_position);
                    obj.error_history(end+1) = error;
                end
                
%             catch ME
%                 warning('����ָ�����ʧ��: %s', ME.message);
            end
        end
        
        % ����SNR
        function snr = estimate_snr(obj, signal)
            try
                if isstruct(signal)
                    if isfield(signal, 'received_signal')
                        signal = signal.received_signal;
                    end
                end
                
                signal = signal(:);
                signal_power = mean(abs(signal).^2);
                
                N = length(signal);
                fft_sig = fft(signal);
                noise_power = mean(abs(fft_sig(floor(N/2):end)).^2);
                
                if noise_power < eps
                    noise_power = eps;
                end
                
                snr = 10*log10(signal_power/noise_power);
                snr = max(min(snr, 60), -20);
                
            catch ME
                fprintf('SNR����ʧ��: %s\n', ME.message);
                snr = obj.sensing_params.min_snr_db;
            end
        end
    end
end