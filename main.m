clear; clc; close all;

start_time = datetime('2025-03-13 11:25:06', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
current_user = 'zz17Pan';
fprintf('程序开始执行时间: %s\n用户: %s\n', datestr(start_time), current_user);

% 初始化系统组件
sensor = MultiModalSensing();
visualizer = SensingVisualizer();
% 初始化动态追踪器（初始状态用当前真实位置及默认角度）
tracker = MultiModelAUKF();

t_sim = 0:0.001:1;  % 1秒仿真，1ms步长

num_frames = length(t_sim);
true_positions = zeros(num_frames, 3);
est_positions = zeros(num_frames, 3);
true_angles = zeros(num_frames, 2);  % [theta, phi]
est_angles = zeros(num_frames, 2);   % [theta, phi]

for idx = 1:num_frames
    t = t_sim(idx);
    fprintf('\n处理第 %d 帧 (t = %.3f s)\n', idx, t);
    % 生成真实轨迹（匀速直线运动）
    true_pos = generate_true_trajectory(t);
    [true_theta, true_phi] = ArrayGeometry.calculate_angles(true_pos);
    % 更新阵列位置，将真实位置传入
    [tx_pos, rx_pos] = ArrayGeometry.initialize_array(t, true_pos);
    % 执行多模态感知
    sensing_result = sensor.perform_multimodal_sensing(tx_pos, rx_pos, t);
    % 更新动态追踪器
    dt = 0.001;
    tracker.predict(dt);
    % 使用测量 [theta; phi; range]（从感知结果获得）
    measurement = [sensing_result.theta; sensing_result.phi; sensing_result.range];
    tracker.update(measurement);
    est_state = tracker.get_estimated_state();
    % 根据追踪器的状态更新估计位置：直接用 x(1:3)
    est_pos = est_state(1:3);
    % 角度取追踪器状态中的 theta, phi
%     est_theta = est_state(7);
%     est_phi = est_state(8);
% 根据估计位置计算角度（替代直接读取状态变量）
    [est_theta, est_phi] = ArrayGeometry.calculate_angles(est_pos);
    

    true_positions(idx,:) = true_pos;
    est_positions(idx,:) = est_pos;
    true_angles(idx,:) = [true_theta, true_phi];
    est_angles(idx,:) = [est_theta, est_phi];
    
    if mod(idx,10)==0
        visualizer.update_plot(true_positions(1:idx,:), est_positions(1:idx,:), true_angles(1:idx,:), est_angles(1:idx,:), t);
        drawnow;
    end
end

visualizer.show_final_results(true_positions, est_positions, true_angles, est_angles, t_sim);
fprintf('\n仿真完成！\n');

%% 生成真实轨迹函数（匀速直线运动）
function pos = generate_true_trajectory(t)
    p0 = [0,0,40];
    v = [1, 0.5, 0];  % 较大运动可放大后观察效果
    pos = p0 + v*t;
end
