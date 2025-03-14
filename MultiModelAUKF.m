classdef MultiModelAUKF < handle
    properties
        x        % 状态向量 [x; y; z; vx; vy; vz]
        P        % 协方差矩阵
        Q        % 过程噪声
        R        % 测量噪声（针对 [theta; phi; range]）
    end
    
    methods
        function obj = MultiModelAUKF(initial_state)
            if nargin < 1
                % 默认初始状态：位置 [0,0,40]，速度 [1,0.5,0]
                initial_state = [0; 0; 40; 1; 0.5; 0];
            end
            obj.x = initial_state;
            obj.P = eye(6)*0.1;
            obj.Q = diag([0.01, 0.01, 0.01, 0.1, 0.1, 0.1]); % 过程噪声
            obj.R = diag([0.001, 0.001, 0.1]); % 测量噪声（角度更敏感）
        end
        
        function predict(obj, dt)
            % CV模型状态转移
            F = eye(6);
            F(1,4) = dt; F(2,5) = dt; F(3,6) = dt;
            obj.x = F * obj.x;
            obj.P = F * obj.P * F' + obj.Q;
        end
        
        function update(obj, measurement)
            % 解析雅可比矩阵（精确推导）
            pos = obj.x(1:3);
            x = pos(1); y = pos(2); z = pos(3);
            r = norm(pos);
            xy_sq = x^2 + y^2;
            
            % 预测测量值
            theta_pred = atan2(y, x);
            phi_pred = atan2(z, sqrt(xy_sq));
            r_pred = r;
            h = [theta_pred; phi_pred; r_pred];
            
            % 雅可比矩阵H
            H = zeros(3,6);
            if abs(xy_sq) > 1e-6
                % d(theta)/dx
                H(1,1) = -y / xy_sq;
                % d(theta)/dy
                H(1,2) = x / xy_sq;
                % d(theta)/dz 及其他列保持0
                
                % d(phi)/dx
                H(2,1) = (-x*z) / (r^2 * sqrt(xy_sq));
                % d(phi)/dy
                H(2,2) = (-y*z) / (r^2 * sqrt(xy_sq));
                % d(phi)/dz
                H(2,3) = sqrt(xy_sq) / r^2;
            end
            
            % d(r)/dx, dy, dz
            H(3,1) = x/r;
            H(3,2) = y/r;
            H(3,3) = z/r;
            
            % 卡尔曼增益计算
            S = H * obj.P * H' + obj.R;
            K = (obj.P * H') / S;
            
            % 状态更新
            innovation = measurement - h;
            obj.x = obj.x + K * innovation;
            obj.P = (eye(6) - K*H) * obj.P;
        end
        
        function state = get_estimated_state(obj)
            state = obj.x;
        end
    end
end