classdef ChannelHSPM < handle
    properties (Constant, Access = private)
        DATE_TIME = '2025-03-13 11:54:21'
        USER = 'zz17Pan'
    end
    
    methods (Static)
        function H = generate_channel(tx_pos, rx_pos, t)
            try
                % 记录开始时间
                start_time = datetime(ChannelHSPM.DATE_TIME, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                fprintf('信道生成开始时间: %s\n用户: %s\n', datestr(start_time), ChannelHSPM.USER);
                
                % 获取配置
                cfg = Config;
                
                % 计算发射接收端总天线数
                [Nt_total, ~] = size(tx_pos);
                [Nr_total, ~] = size(rx_pos);
                
                fprintf('初始化信道矩阵，维度: [%d, %d]\n', Nr_total, Nt_total);
                
                % 初始化HSPM信道矩阵
                H = zeros(Nr_total, Nt_total);
                
                % 子阵天线数
                N_at = cfg.Nx * cfg.Nz;
                
                % 对每个子阵对生成信道矩阵
                for kr = 1:cfg.Kr
                    for kt = 1:cfg.Kt
                        % 获取当前子阵的天线索引
                        tx_indices = ((kt-1)*N_at + 1):(kt*N_at);
                        rx_indices = ((kr-1)*N_at + 1):(kr*N_at);
                        
                        % 获取子阵天线位置
                        tx_block = tx_pos(tx_indices, :);
                        rx_block = rx_pos(rx_indices, :);
                        
                        % 计算子阵中心
                        tx_center = mean(tx_block);
                        rx_center = mean(rx_block);
                        
                        % 计算平均距离
                        D_mean = norm(rx_center - tx_center);
                        
                        % 计算LoS因子
                        LoS_factor = ChannelHSPM.calculate_los_factor(D_mean, cfg);
                        
                        % 计算角度
                        [theta_rx, phi_rx] = ChannelHSPM.calculate_angles(rx_center - tx_center);
                        [theta_tx, phi_tx] = ChannelHSPM.calculate_angles(tx_center - rx_center);
                        
                        % 生成子阵方向矢量
                        a_rx = ChannelHSPM.steering_vector_2D(theta_rx, phi_rx, cfg.Nx, cfg.Nz, cfg.lambda, cfg.d);
                        a_tx = ChannelHSPM.steering_vector_2D(theta_tx, phi_tx, cfg.Nx, cfg.Nz, cfg.lambda, cfg.d);
                        
                        % 生成极化矩阵
                        P = ChannelHSPM.generate_polarization_matrix(theta_rx, phi_rx, theta_tx, phi_tx);
                        
                        % 生成子阵信道块（利用 Kronecker 扩展极化矩阵）
                        H_pol = kron(P, ones(N_at/2, N_at/2));
                        H_block = LoS_factor * (a_rx * a_tx') .* H_pol;
                        
                        % 填充到总信道矩阵
                        H(rx_indices, tx_indices) = H_block;
                        
                        fprintf('处理子阵 [%d,%d]: 距离=%.2fm, θ_rx=%.2f°, φ_rx=%.2f°\n', ...
                            kr, kt, D_mean, rad2deg(theta_rx), rad2deg(phi_rx));
                    end
                end
                
                % 应用功率缩放
                H = cfg.Beta * H;
                
                % 添加噪声
                noise_power = 10^(-cfg.SNR/10);
                H = H + sqrt(noise_power/2) * (randn(size(H)) + 1j*randn(size(H)));
                
                % 验证信道矩阵
                ChannelHSPM.verify_channel_matrix(H);
                
                fprintf('信道生成完成\n');
                
            catch ME
                fprintf('信道生成失败：\n');
                fprintf('错误信息：%s\n', ME.message);
                fprintf('错误位置：%s\n', ME.stack(1).name);
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
                fprintf('方向矢量生成失败：\n');
                fprintf('错误信息：%s\n', ME.message);
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
                fprintf('极化矩阵生成失败：\n');
                fprintf('错误信息：%s\n', ME.message);
                rethrow(ME);
            end
        end
        
        function verify_channel_matrix(H)
            try
                if any(size(H) <= 0)
                    error('信道矩阵维度无效');
                end
                if any(isnan(H(:))) || any(isinf(H(:)))
                    error('信道矩阵包含无效值');
                end
                cond_num = cond(H);
                if cond_num > 1e6
                    warning('信道矩阵条件数过大: %.2e', cond_num);
                end
                fprintf('信道矩阵统计：\n');
                fprintf('维度: [%d, %d]\n', size(H));
                fprintf('均值: %.2e + %.2ej\n', mean(real(H(:))), mean(imag(H(:))));
                fprintf('标准差: %.2e\n', std(abs(H(:))));
                fprintf('条件数: %.2e\n', cond_num);
            catch ME
                fprintf('信道矩阵验证失败：\n');
                fprintf('错误信息：%s\n', ME.message);
                rethrow(ME);
            end
        end
    end
end
