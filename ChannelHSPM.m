classdef ChannelHSPM < handle
    properties (Constant, Access = private)
        DATE_TIME = '2025-03-13 11:54:21'
        USER = 'zz17Pan'
    end
    
    methods (Static)
        function H = generate_channel(tx_pos, rx_pos, t)
            try
                % ��¼��ʼʱ��
                start_time = datetime(ChannelHSPM.DATE_TIME, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                fprintf('�ŵ����ɿ�ʼʱ��: %s\n�û�: %s\n', datestr(start_time), ChannelHSPM.USER);
                
                % ��ȡ����
                cfg = Config;
                
                % ���㷢����ն���������
                [Nt_total, ~] = size(tx_pos);
                [Nr_total, ~] = size(rx_pos);
                
                fprintf('��ʼ���ŵ�����ά��: [%d, %d]\n', Nr_total, Nt_total);
                
                % ��ʼ��HSPM�ŵ�����
                H = zeros(Nr_total, Nt_total);
                
                % ����������
                N_at = cfg.Nx * cfg.Nz;
                
                % ��ÿ������������ŵ�����
                for kr = 1:cfg.Kr
                    for kt = 1:cfg.Kt
                        % ��ȡ��ǰ�������������
                        tx_indices = ((kt-1)*N_at + 1):(kt*N_at);
                        rx_indices = ((kr-1)*N_at + 1):(kr*N_at);
                        
                        % ��ȡ��������λ��
                        tx_block = tx_pos(tx_indices, :);
                        rx_block = rx_pos(rx_indices, :);
                        
                        % ������������
                        tx_center = mean(tx_block);
                        rx_center = mean(rx_block);
                        
                        % ����ƽ������
                        D_mean = norm(rx_center - tx_center);
                        
                        % ����LoS����
                        LoS_factor = ChannelHSPM.calculate_los_factor(D_mean, cfg);
                        
                        % ����Ƕ�
                        [theta_rx, phi_rx] = ChannelHSPM.calculate_angles(rx_center - tx_center);
                        [theta_tx, phi_tx] = ChannelHSPM.calculate_angles(tx_center - rx_center);
                        
                        % ����������ʸ��
                        a_rx = ChannelHSPM.steering_vector_2D(theta_rx, phi_rx, cfg.Nx, cfg.Nz, cfg.lambda, cfg.d);
                        a_tx = ChannelHSPM.steering_vector_2D(theta_tx, phi_tx, cfg.Nx, cfg.Nz, cfg.lambda, cfg.d);
                        
                        % ���ɼ�������
                        P = ChannelHSPM.generate_polarization_matrix(theta_rx, phi_rx, theta_tx, phi_tx);
                        
                        % ���������ŵ��飨���� Kronecker ��չ��������
                        H_pol = kron(P, ones(N_at/2, N_at/2));
                        H_block = LoS_factor * (a_rx * a_tx') .* H_pol;
                        
                        % ��䵽���ŵ�����
                        H(rx_indices, tx_indices) = H_block;
                        
                        fprintf('�������� [%d,%d]: ����=%.2fm, ��_rx=%.2f��, ��_rx=%.2f��\n', ...
                            kr, kt, D_mean, rad2deg(theta_rx), rad2deg(phi_rx));
                    end
                end
                
                % Ӧ�ù�������
                H = cfg.Beta * H;
                
                % �������
                noise_power = 10^(-cfg.SNR/10);
                H = H + sqrt(noise_power/2) * (randn(size(H)) + 1j*randn(size(H)));
                
                % ��֤�ŵ�����
                ChannelHSPM.verify_channel_matrix(H);
                
                fprintf('�ŵ��������\n');
                
            catch ME
                fprintf('�ŵ�����ʧ�ܣ�\n');
                fprintf('������Ϣ��%s\n', ME.message);
                fprintf('����λ�ã�%s\n', ME.stack(1).name);
                rethrow(ME);
            end
        end
        
        function [theta, phi] = calculate_angles(pos)
            x = pos(1); y = pos(2); z = pos(3);
            r_xy = sqrt(x^2 + y^2);
            theta = atan2(y, x);
            phi = atan2(z, r_xy);
            
            if r_xy < 1e-10
                theta = 0;
                phi = sign(z) * pi/2;
            end
            
            theta = mod(theta + pi, 2*pi) - pi;
            phi = max(min(phi, pi/2), -pi/2);
        end
        
        function a = steering_vector_2D(theta, phi, Nx, Nz, lambda, d)
            try
                a = zeros(Nx*Nz, 1);
                k = 2*pi/lambda;
                idx = 1;
                for nz = 0:(Nz-1)
                    for nx = 0:(Nx-1)
                        phase = k * d * (nx * sin(theta) * cos(phi) + nz * sin(phi));
                        a(idx) = exp(1j*phase);
                        idx = idx + 1;
                    end
                end
                a = a / norm(a);
            catch ME
                fprintf('����ʸ������ʧ�ܣ�\n');
                fprintf('������Ϣ��%s\n', ME.message);
                rethrow(ME);
            end
        end
        
        function los_factor = calculate_los_factor(D, cfg)
            path_loss = cfg.c / (4*pi*cfg.f*D);
            attenuation = exp(-0.5*cfg.k_f*D);
            phase_rotation = exp(-1j*2*pi*D/cfg.lambda);
            los_factor = path_loss * attenuation * phase_rotation;
        end
        
        function P = generate_polarization_matrix(theta_rx, phi_rx, theta_tx, phi_tx)
            try
                e_theta_rx = [-sin(theta_rx); cos(theta_rx); 0];
                e_phi_rx = [-cos(theta_rx)*sin(phi_rx); -sin(theta_rx)*sin(phi_rx); cos(phi_rx)];
                e_theta_tx = [-sin(theta_tx); cos(theta_tx); 0];
                e_phi_tx = [-cos(theta_tx)*sin(phi_tx); -sin(theta_tx)*sin(phi_tx); cos(phi_tx)];
                P = [dot(e_theta_rx, e_theta_tx), dot(e_theta_rx, e_phi_tx);
                     dot(e_phi_rx, e_theta_tx), dot(e_phi_rx, e_phi_tx)];
                P = P / norm(P, 'fro');
            catch ME
                fprintf('������������ʧ�ܣ�\n');
                fprintf('������Ϣ��%s\n', ME.message);
                rethrow(ME);
            end
        end
        
        function verify_channel_matrix(H)
            try
                if any(size(H) <= 0)
                    error('�ŵ�����ά����Ч');
                end
                if any(isnan(H(:))) || any(isinf(H(:)))
                    error('�ŵ����������Чֵ');
                end
                cond_num = cond(H);
                if cond_num > 1e6
                    warning('�ŵ���������������: %.2e', cond_num);
                end
                fprintf('�ŵ�����ͳ�ƣ�\n');
                fprintf('ά��: [%d, %d]\n', size(H));
                fprintf('��ֵ: %.2e + %.2ej\n', mean(real(H(:))), mean(imag(H(:))));
                fprintf('��׼��: %.2e\n', std(abs(H(:))));
                fprintf('������: %.2e\n', cond_num);
            catch ME
                fprintf('�ŵ�������֤ʧ�ܣ�\n');
                fprintf('������Ϣ��%s\n', ME.message);
                rethrow(ME);
            end
        end
    end
end
