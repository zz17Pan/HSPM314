classdef ArrayGeometry < handle 
    methods (Static)
        % ��������ʵλ�� true_rx ���ڸ��� rx_centers
        function [tx_pos, rx_pos] = initialize_array(t, true_rx)
            try
                cfg = Config;
                [tx_centers, rx_centers] = ArrayGeometry.calculate_subarray_centers(t, true_rx);
                
                % ���ɷ����λ�þ��󣨹̶��ṹ��
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
                
                % ���ɽ��ն�λ�þ����ڹ̶�����ṹ�����ϼ�����ʵ������λ��
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
                
                fprintf('���г�ʼ����� (t = %.3f s):\n', t);
                fprintf('��������ά��: [%d, %d]\n', size(tx_pos));
                fprintf('��������ά��: [%d, %d]\n', size(rx_pos));
            catch ME
                fprintf('���г�ʼ��ʧ�ܣ�\n%s\n', ME.message);
                rethrow(ME);
            end
        end
        
        % ����ʵ������λ�ü������
        function [tx_centers, rx_centers] = calculate_subarray_centers(t, true_rx)
            cfg = Config;
            % ����˹̶���ԭ��
            tx_centers = zeros(cfg.Kt, 3);
            for kt = 1:cfg.Kt
                tx_centers(kt,1) = (kt - (cfg.Kt+1)/2) * cfg.d_sub;
            end
            % ���ն����������ڹ̶����ֻ����ϼ�����ʵ������λ��
            rx_centers = zeros(cfg.Kr, 3);
            for kr = 1:cfg.Kr
                % ���������� x �����У�����ڽ��������ĵ�ƫ��
                rx_center_fixed = [(kr - (cfg.Kr+1)/2) * cfg.d_sub, 0, 0];
                rx_centers(kr,:) = true_rx + rx_center_fixed;
            end
        end
        
        function [theta, phi] = calculate_angles(pos)
            x = pos(1); y = pos(2); z = pos(3);
            r_xy = sqrt(x^2 + y^2);
            theta = atan2(y, x);
            phi = atan2(z, r_xy);
            fprintf('�Ƕȼ�����: ��λ��=%.2f�� (%.4f rad), ������=%.2f�� (%.4f rad)\n', rad2deg(theta), theta, rad2deg(phi), phi);
        end
    end
end
