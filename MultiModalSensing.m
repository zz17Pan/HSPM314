classdef MultiModalSensing < handle
    properties
        cfg                 % 系统配置参数
        dict               % 联合字典
        sensing_params     % 感知参数
        sensor_state      % 当前感知状态
        valid_frames = 0   % 有效帧计数
        snr_history = []   % SNR历史记录
        error_history = [] % 误差历史记录
    end
    
    methods
        % 构造函数
        function obj = MultiModalSensing()
            try
                obj.cfg = Config;
                
                % 设置感知参数
                obj.sensing_params = struct(...
                    'fmcw_bandwidth', 1e9, ...    % 1GHz带宽
                    'sweep_time', 1e-3, ...       % 1ms扫描时间
                    'N_samples', 1024, ...        % 采样点数
                    'N_range_fft', 2048, ...      % 距离FFT点数
                    'N_doppler_fft', 256, ...     % 多普勒FFT点数
                    'dict_size', [180, 180, 50],... % 字典维度
                    'omp_sparsity', 3, ...        % OMP稀疏度
                    'min_snr_db', 5);             % 最小SNR阈值
                
                fprintf('系统参数: Nx=%d, Nz=%d, 天线总数=%d\n', ...
                    obj.cfg.Nx, obj.cfg.Nz, obj.cfg.Nx * obj.cfg.Nz);
                
                fprintf('开始生成联合字典...\n');
                obj.dict = obj.generate_dictionary();
                fprintf('字典生成完成\n');
                
            catch ME
                fprintf('初始化失败: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % 多模态感知主函数
        function result = perform_multimodal_sensing(obj, tx_pos, rx_pos, t)
            try
                fprintf('\n处理第 %d 帧 (t = %.3f s)\n', obj.valid_frames + 1, t);
                
                % 1. 角度域感知
                angle_result = obj.perform_angle_sensing(tx_pos, rx_pos, t);
                
                % 验证结果完整性
                required_fields = {'position', 'theta', 'phi', 'range', 'received_signal'};
                for i = 1:length(required_fields)
                    if ~isfield(angle_result, required_fields{i})
                        error('角度感知结果缺少必要字段: %s', required_fields{i});
                    end
                end
                
                % 2. 距离-多普勒感知
                range_doppler = obj.perform_range_doppler_sensing(tx_pos, rx_pos);
                
                % 3. 极化域感知
                polarization = obj.perform_polarization_sensing(tx_pos, rx_pos);
                
                % 4. 特征融合
                fusion_result = obj.feature_fusion(angle_result, range_doppler, polarization);
                
                % 构造完整结果
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
                fprintf('多模态感知处理失败: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % 角度域感知
        function result = perform_angle_sensing(obj, tx_pos, rx_pos, t)
            try
                % 1. 获取子阵配置
                [~, sensing_subarrays] = ResourceAllocation.allocate_subarrays();
                
                % 2. 生成导频信号
                N_tx = size(tx_pos, 1);
                pilot = obj.generate_pilot_signal(N_tx);
                
                % 3. 生成信道矩阵
                H = ChannelHSPM.generate_channel(tx_pos, rx_pos, t);
                
                % 4. 接收信号处理
                Y = H * pilot;
                
                % 5. 维度匹配处理
                N_ant = obj.cfg.Nx * obj.cfg.Nz;
                if size(Y,1) > N_ant
                    Y = Y(1:N_ant, :);
                end
                
                % 6. SNR估计
                snr = obj.estimate_snr(Y);
                fprintf('估计SNR: %.2f dB\n', snr);
                
                % 7. 参数估计
                if snr >= obj.sensing_params.min_snr_db
                    [theta, phi, r] = obj.estimate_parameters(Y, obj.dict);
                else
                    fprintf('SNR低于阈值，使用上一帧结果\n');
                    last_est = obj.get_last_estimate();
                    theta = last_est.theta;
                    phi = last_est.phi;
                    r = last_est.range;
                end
                
                % 8. 构造结果
                result = struct(...
                    'position', [r*cos(phi)*cos(theta), r*cos(phi)*sin(theta), r*sin(phi)], ...
                    'theta', theta, ...
                    'phi', phi, ...
                    'range', r, ...
                    'velocity', 0, ...
                    'received_signal', Y);
                
                fprintf('角度估计结果：\n');
                fprintf('方位角: %.2f°\n', rad2deg(theta));
                fprintf('俯仰角: %.2f°\n', rad2deg(phi));
                fprintf('距离: %.2f m\n', r);
                
            catch ME
                fprintf('角度感知失败: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % 距离多普勒感知
        function result = perform_range_doppler_sensing(obj, tx_pos, rx_pos)
            try
                % 1. FMCW参数设置
                params = obj.sensing_params;
                N_samples = params.N_samples;
                Ts = params.sweep_time / N_samples;
                t_samples = (0:N_samples-1) * Ts;
                
                % 2. 生成FMCW信号
                fmcw = obj.generate_fmcw_signal(t_samples, params.fmcw_bandwidth, params.sweep_time);
                
                % 3. 生成信道矩阵
                H = ChannelHSPM.generate_channel(tx_pos, rx_pos, 0);
                
                % 4. 接收信号处理
                tx_signal = repmat(fmcw(:), 1, size(H, 2))';
                rx_signal = H * tx_signal;
                
                % 5. 初始化距离和多普勒剖面
                N_rx = size(rx_signal, 1);
                range_profiles = zeros(N_rx, params.N_range_fft/2);
                doppler_profiles = zeros(N_rx, params.N_doppler_fft/2);
                
                % 6. 处理每个天线信号
                for n = 1:N_rx
                    [range_prof, doppler_prof] = obj.process_fmcw(...
                        rx_signal(n,:), params.N_range_fft, params.N_doppler_fft);
                    range_profiles(n,:) = range_prof;
                    doppler_profiles(n,:) = doppler_prof;
                end
                
                % 7. 构造结果
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
                fprintf('距离多普勒感知失败: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % 生成FMCW信号
        function signal = generate_fmcw_signal(~, t, bandwidth, sweep_time)
            % 计算调频斜率
            mu = bandwidth / sweep_time;
            
            % 生成基带FMCW信号
            signal = exp(1j*pi*mu*t.^2);
            
            % 添加窗函数和归一化
            window = hamming(length(t));
            signal = signal(:) .* window;
            signal = signal / norm(signal);
        end

        % 极化域感知
        function result = perform_polarization_sensing(obj, tx_pos, rx_pos)
            try
                % 获取极化信道
                H = ChannelHSPM.generate_channel(tx_pos, rx_pos, 0);
                
                % 提取极化分量
                N_ant = obj.cfg.Nx * obj.cfg.Nz;
                H_HH = H(1:N_ant/2, 1:N_ant/2);
                H_HV = H(1:N_ant/2, N_ant/2+1:end);
                H_VH = H(N_ant/2+1:end, 1:N_ant/2);
                H_VV = H(N_ant/2+1:end, N_ant/2+1:end);
                
                % 构造极化矩阵
                H_pol = [mean2(H_HH) mean2(H_HV); mean2(H_VH) mean2(H_VV)];
                
                % 估计极化参数
                [U, S, V] = svd(H_pol);
                pol_angle = angle(V(1,1));
                ellipticity = S(2,2) / S(1,1);
                
                result = struct(...
                    'pol_angle', pol_angle, ...
                    'ellipticity', ellipticity, ...
                    'received_signal', H);
                
            catch ME
                fprintf('极化感知失败: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % 特征融合
        function result = feature_fusion(obj, angle_result, range_doppler, polarization)
            try
                % 计算各测量置信度
                conf_angle = obj.calculate_confidence(obj.estimate_snr(angle_result.received_signal));
                conf_range = obj.calculate_confidence(obj.estimate_snr(range_doppler.received_signal));
                conf_pol = obj.calculate_confidence(obj.estimate_snr(polarization.received_signal));
                
                % D-S证据理论融合
                weights = [conf_angle, conf_range, conf_pol];
                fusion_coef = obj.ds_fusion(weights);
                
                % 位置融合
                pos_angle = angle_result.position;
                pos_range = [range_doppler.range * cos(angle_result.theta), ...
                           range_doppler.range * sin(angle_result.theta), ...
                           range_doppler.range * sin(angle_result.phi)];
                           
                fused_position = fusion_coef * pos_angle + (1-fusion_coef) * pos_range;
                
                result = struct(...
                    'position', fused_position, ...
                    'confidence', fusion_coef);
                
            catch ME
                fprintf('特征融合失败: %s\n', ME.message);
                rethrow(ME);
            end
        end

        % 生成字典
        function dict = generate_dictionary(obj)
            try
                % 1. 获取系统参数
                Nx = obj.cfg.Nx;
                Nz = obj.cfg.Nz;
                N_ant = Nx * Nz;
                dict_size = obj.sensing_params.dict_size;
                total_atoms = dict_size(1) * dict_size(2) * dict_size(3);
                
                % 2. 初始化字典矩阵
                dict = zeros(N_ant, total_atoms);
                fprintf('开始生成字典，原子总数: %d\n', total_atoms);
                
                % 3. 生成参数网格
                theta_grid = linspace(-pi/2, pi/2, dict_size(1));
                phi_grid = linspace(-pi/2, pi/2, dict_size(2));
                r_grid = linspace(1, 100, dict_size(3));
                
                % 4. 生成字典原子
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
                                fprintf('字典生成进度: %.1f%%\n', 100*idx/total_atoms);
                            end
                            idx = idx + 1;
                        end
                    end
                end
                
                % 5. 正交化和归一化
                fprintf('正交化处理...\n');
                [dict, ~] = qr(dict, 0);
                for i = 1:size(dict, 2)
                    dict(:,i) = dict(:,i) / norm(dict(:,i));
                end
                
                fprintf('字典维度: [%d, %d]\n', size(dict));
                
            catch ME
                fprintf('字典生成失败: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
%         % 生成字典原子
%         function atom = generate_atom(obj, theta, phi, r)
%             try
%                 % 1. 获取系统参数
%                 Nx = obj.cfg.Nx;
%                 Nz = obj.cfg.Nz;
%                 d = obj.cfg.d;
%                 lambda = obj.cfg.lambda;
%                 k = 2*pi/lambda;
%                                 % 继续 generate_atom 函数
%                 % 2. 初始化原子向量
%                 atom = zeros(Nx*Nz, 1);
%                 idx = 1;
%                 
%                 % 3. 生成空间响应
%                 for nz = 0:(Nz-1)
%                     for nx = 0:(Nx-1)
%                         % 计算相位
%                         phase = k * d * (...
%                             nx * sin(theta) * cos(phi) + ...
%                             nz * sin(phi));
%                         
%                         % 添加距离衰减
%                         amp = exp(-r/100);
%                         
%                         % 生成复响应
%                         atom(idx) = amp * exp(1j*phase);
%                         idx = idx + 1;
%                     end
%                 end
%                 
%                 % 4. 归一化
%                 atom = atom / norm(atom);
%                 
%             catch ME
%                 fprintf('原子生成失败: %s\n', ME.message);
%                 rethrow(ME);
%             end
%         end
function atom = generate_atom(obj, theta, phi, r)
    % 1. 获取系统参数
    Nx = obj.cfg.Nx;
    Nz = obj.cfg.Nz;
    d = obj.cfg.d;
    lambda = obj.cfg.lambda;
    k = 2*pi/lambda;
    
    % 2. 自由空间路径损耗模型
    path_loss = (lambda/(4*pi*r))^2;  % 功率衰减
    amp = sqrt(path_loss);            % 幅度衰减
    
    % 3. 初始化原子向量
    atom = zeros(Nx*Nz, 1);
    idx = 1;
    
    % 4. 生成空间响应
    for nz = 0:(Nz-1)
        for nx = 0:(Nx-1)
            % 相位计算（考虑三维几何）
            phase = k * d * (...
                nx * sin(theta) * cos(phi) + ...
                nz * sin(phi)...
            );
            atom(idx) = amp * exp(1j*phase);
            idx = idx + 1;
        end
    end
    
    % 5. 归一化
    atom = atom / norm(atom);
end
        % OMP算法实现
        function x = omp(obj, dict, y, K)
            try
                % 1. 检查和预处理输入
                [M, N] = size(dict);
                y = y(:);  % 确保y是列向量
                
                if length(y) ~= M
                    error('维度不匹配: dict(%d×%d), y(%d×1)', M, N, length(y));
                end
                
                % 2. 初始化
                x = zeros(N, 1);
                residual = y;
                support = [];
                D_s = [];
                x_s = [];
                
                % 3. 迭代重建
                for k = 1:min(K, N)
                    % 计算相关性
                    corr = abs(dict' * residual);
                    [~, idx] = max(corr);
                    
                    % 更新支持集
                    if ~ismember(idx, support)
                        support = [support; idx];
                        
                        % 最小二乘估计
                        D_s = dict(:, support);
                        if rank(D_s) < length(support)
                            support = support(1:end-1);
                            continue;
                        end
                        
                        x_s = pinv(D_s) * y;
                        
                        % 更新残差
                        residual = y - D_s * x_s;
                        
                        % 收敛检查
                        if norm(residual) < 1e-6
                            break;
                        end
                    end
                end
                
                % 4. 填充结果
                x(support) = x_s;
                
            catch ME
                fprintf('OMP算法失败: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % 参数估计
        function [theta, phi, r] = estimate_parameters(obj, Y, dict)
            try
                % 1. 检查输入维度
                [M1, N1] = size(dict);
                [M2, N2] = size(Y);
                
                % 2. 确保Y是列向量
                if N2 > 1
                    Y = Y(:,1);
                end
                
                % 3. 验证维度匹配
                if M1 ~= M2
                    error('字典和观测信号维度不匹配：dict(%d×%d), Y(%d×%d)', ...
                        M1, N1, M2, N2);
                end
                
                % 4. 执行OMP重建
                x = obj.omp(dict, Y, obj.sensing_params.omp_sparsity);
                
                % 5. 找到最强分量
                [~, max_idx] = max(abs(x));
                
                % 6. 从索引恢复参数
                dict_size = obj.sensing_params.dict_size;
                [i_theta, i_phi, i_r] = ind2sub(dict_size, max_idx);
                
                % 7. 计算实际参数值
                theta_grid = linspace(-pi/2, pi/2, dict_size(1));
                phi_grid = linspace(-pi/2, pi/2, dict_size(2));
                r_grid = linspace(1, 100, dict_size(3));
                
                theta = theta_grid(i_theta);
                phi = phi_grid(i_phi);
                r = r_grid(i_r);
                
            catch ME
                fprintf('参数估计失败: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % 生成导频信号
        function pilot = generate_pilot_signal(~, N_tx)
            try
                pilot = eye(N_tx);
                phases = exp(1j*2*pi*rand(N_tx, 1));
                pilot = pilot * diag(phases);
                pilot = pilot / sqrt(N_tx);
                
            catch ME
                fprintf('导频信号生成失败: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % FMCW信号处理
        function [range_profile, doppler_profile] = process_fmcw(~, rx_signal, N_range_fft, N_doppler_fft)
            try
                % 1. 确保信号为行向量
                rx_signal = rx_signal(:).';
                
                % 2. 应用窗函数
                window = hamming(length(rx_signal));
                rx_signal_windowed = rx_signal .* window.';
                
                % 3. 距离FFT处理
                range_fft = fft(rx_signal_windowed, N_range_fft);
                range_profile = abs(range_fft(1:N_range_fft/2));
                
                % 4. 多普勒FFT处理
                doppler_fft = fft(rx_signal_windowed, N_doppler_fft);
                doppler_profile = abs(doppler_fft(1:N_doppler_fft/2));
                
                % 5. 归一化
                range_profile = range_profile / max(abs(range_profile) + eps);
                doppler_profile = doppler_profile / max(abs(doppler_profile) + eps);
                
            catch ME
                fprintf('FMCW处理失败: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % 计算置信度
        function conf = calculate_confidence(~, SNR)
            conf = 1 - exp(-SNR/10);
        end
        
        % D-S证据理论融合
        function result = ds_fusion(~, masses)
            result = prod(masses) / (prod(masses) + prod(1-masses));
        end
        
        % 获取上一帧估计结果
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
                fprintf('获取上一帧估计失败，使用默认值\n');
                last_est = struct(...
                    'position', [0, 0, 40], ...
                    'theta', 0, ...
                    'phi', pi/4, ...
                    'range', 40, ...
                    'velocity', 0);
            end
        end
        
        % 更新性能指标
        function update_performance_metrics(obj, result)
            try
                current_snr = obj.estimate_snr(result.received_signal);
                obj.snr_history(end+1) = current_snr;
                
                if isfield(result, 'true_position') && isfield(result, 'position')
                    error = norm(result.position - result.true_position);
                    obj.error_history(end+1) = error;
                end
                
%             catch ME
%                 warning('性能指标更新失败: %s', ME.message);
            end
        end
        
        % 估计SNR
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
                fprintf('SNR估计失败: %s\n', ME.message);
                snr = obj.sensing_params.min_snr_db;
            end
        end
    end
end