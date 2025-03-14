classdef ArrayGeometry < handle 
    methods (Static)
        % 接收器真实位置 true_rx 用于更新 rx_centers
        function [tx_pos, rx_pos] = initialize_array(t, true_rx)
            try
                cfg = Config;
                [tx_centers, rx_centers] = ArrayGeometry.calculate_subarray_centers(t, true_rx);
                
                % 生成发射端位置矩阵（固定结构）
                tx_pos = zeros(cfg.Kt * cfg.Nx * cfg.Nz, 3);
                idx_tx = 1;
                for kt = 1:cfg.Kt
                    for nz = 0:(cfg.Nz-1)
                        for nx = 0:(cfg.Nx-1)
                            local_x = (nx - (cfg.Nx-1)/2) * cfg.d;
                            local_z = (nz - (cfg.Nz-1)/2) * cfg.d;
                            tx_pos(idx_tx,:) = tx_centers(kt,:) + [local_x, 0, local_z];
                            idx_tx = idx_tx + 1;
                        end
                    end
                end
                
                % 生成接收端位置矩阵：在固定子阵结构基础上加上真实接收器位置
                rx_pos = zeros(cfg.Kr * cfg.Nx * cfg.Nz, 3);
                idx_rx = 1;
                for kr = 1:cfg.Kr
                    for nz = 0:(cfg.Nz-1)
                        for nx = 0:(cfg.Nx-1)
                            local_x = (nx - (cfg.Nx-1)/2) * cfg.d;
                            local_z = (nz - (cfg.Nz-1)/2) * cfg.d;
                            rx_pos(idx_rx,:) = rx_centers(kr,:) + [local_x, 0, local_z];
                            idx_rx = idx_rx + 1;
                        end
                    end
                end
                
                fprintf('阵列初始化完成 (t = %.3f s):\n', t);
                fprintf('发射阵列维度: [%d, %d]\n', size(tx_pos));
                fprintf('接收阵列维度: [%d, %d]\n', size(rx_pos));
            catch ME
                fprintf('阵列初始化失败：\n%s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % 将真实接收器位置加入计算
        function [tx_centers, rx_centers] = calculate_subarray_centers(t, true_rx)
            cfg = Config;
            % 发射端固定于原点
            tx_centers = zeros(cfg.Kt, 3);
            for kt = 1:cfg.Kt
                tx_centers(kt,1) = (kt - (cfg.Kt+1)/2) * cfg.d_sub;
            end
            % 接收端子阵中心在固定布局基础上加上真实接收器位置
            rx_centers = zeros(cfg.Kr, 3);
            for kr = 1:cfg.Kr
                % 假设子阵沿 x 轴排列，相对于接收器中心的偏移
                rx_center_fixed = [(kr - (cfg.Kr+1)/2) * cfg.d_sub, 0, 0];
                rx_centers(kr,:) = true_rx + rx_center_fixed;
            end
        end
        
        function [theta, phi] = calculate_angles(pos)
            x = pos(1); y = pos(2); z = pos(3);
            r_xy = sqrt(x^2 + y^2);
            theta = atan2(y, x);
            phi = atan2(z, r_xy);
            fprintf('角度计算结果: 方位角=%.2f° (%.4f rad), 俯仰角=%.2f° (%.4f rad)\n', rad2deg(theta), theta, rad2deg(phi), phi);
        end
    end
end
