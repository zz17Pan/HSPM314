clear; clc; close all;

start_time = datetime('2025-03-13 11:25:06', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
current_user = 'zz17Pan';
fprintf('����ʼִ��ʱ��: %s\n�û�: %s\n', datestr(start_time), current_user);

% ��ʼ��ϵͳ���
sensor = MultiModalSensing();
visualizer = SensingVisualizer();
% ��ʼ����̬׷��������ʼ״̬�õ�ǰ��ʵλ�ü�Ĭ�ϽǶȣ�
tracker = MultiModelAUKF();

t_sim = 0:0.001:1;  % 1����棬1ms����

num_frames = length(t_sim);
true_positions = zeros(num_frames, 3);
est_positions = zeros(num_frames, 3);
true_angles = zeros(num_frames, 2);  % [theta, phi]
est_angles = zeros(num_frames, 2);   % [theta, phi]

for idx = 1:num_frames
    t = t_sim(idx);
    fprintf('\n����� %d ֡ (t = %.3f s)\n', idx, t);
    % ������ʵ�켣������ֱ���˶���
    true_pos = generate_true_trajectory(t);
    [true_theta, true_phi] = ArrayGeometry.calculate_angles(true_pos);
    % ��������λ�ã�����ʵλ�ô���
    [tx_pos, rx_pos] = ArrayGeometry.initialize_array(t, true_pos);
    % ִ�ж�ģ̬��֪
    sensing_result = sensor.perform_multimodal_sensing(tx_pos, rx_pos, t);
    % ���¶�̬׷����
    dt = 0.001;
    tracker.predict(dt);
    % ʹ�ò��� [theta; phi; range]���Ӹ�֪�����ã�
    measurement = [sensing_result.theta; sensing_result.phi; sensing_result.range];
    tracker.update(measurement);
    est_state = tracker.get_estimated_state();
    % ����׷������״̬���¹���λ�ã�ֱ���� x(1:3)
    est_pos = est_state(1:3);
    % �Ƕ�ȡ׷����״̬�е� theta, phi
%     est_theta = est_state(7);
%     est_phi = est_state(8);
% ���ݹ���λ�ü���Ƕȣ����ֱ�Ӷ�ȡ״̬������
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
fprintf('\n������ɣ�\n');

%% ������ʵ�켣����������ֱ���˶���
function pos = generate_true_trajectory(t)
    p0 = [0,0,40];
    v = [1, 0.5, 0];  % �ϴ��˶��ɷŴ��۲�Ч��
    pos = p0 + v*t;
end
