classdef SensingVisualizer < handle
    properties
        fig_handle
        subplot_handles
    end
    
    methods
        function obj = SensingVisualizer()
            obj.fig_handle = figure('Name', '��ģ̬��֪ʵʱ���', 'Position', [100 100 1200 800]);
            obj.subplot_handles.position = subplot(2,2,1);
            obj.subplot_handles.angles = subplot(2,2,2);
            obj.subplot_handles.error = subplot(2,2,3);
            obj.subplot_handles.stats = subplot(2,2,4);
        end
        
        function update_plot(obj, true_pos, est_pos, true_ang, est_ang, t)
            subplot(obj.subplot_handles.position);
            plot3(true_pos(:,1), true_pos(:,2), true_pos(:,3), 'b-', 'LineWidth', 2);
            hold on;
            plot3(est_pos(:,1), est_pos(:,2), est_pos(:,3), 'r--', 'LineWidth', 2);
            scatter3(true_pos(end,1), true_pos(end,2), true_pos(end,3), 100, 'b', 'filled');
            scatter3(est_pos(end,1), est_pos(end,2), est_pos(end,3), 100, 'r', 'filled');
            hold off;
            grid on;
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
            legend('��ʵ�켣', '���ƹ켣', '��ʵλ��', '����λ��');
            title(sprintf('3Dλ�ø��� (t = %.3fs)', t));
            
            subplot(obj.subplot_handles.angles);
            plot(true_ang(:,1)*180/pi, true_ang(:,2)*180/pi, 'b-', 'LineWidth', 2);
            hold on;
            plot(est_ang(:,1)*180/pi, est_ang(:,2)*180/pi, 'r--', 'LineWidth', 2);
            scatter(true_ang(end,1)*180/pi, true_ang(end,2)*180/pi, 100, 'b', 'filled');
            scatter(est_ang(end,1)*180/pi, est_ang(end,2)*180/pi, 100, 'r', 'filled');
            hold off;
            grid on;
            xlabel('��λ�� (��)'); ylabel('������ (��)');
            legend('��ʵ�Ƕ�', '���ƽǶ�', '��ʵ��ǰ', '���Ƶ�ǰ');
            title('�Ƕȸ���');
            
            subplot(obj.subplot_handles.error);
            pos_error = sqrt(sum((true_pos - est_pos).^2,2));
            ang_error = sqrt(sum((true_ang - est_ang).^2,2));
            plot(pos_error, 'b-', 'LineWidth', 1.5);
            hold on;
            plot(ang_error*180/pi, 'r--', 'LineWidth', 1.5);
            hold off;
            grid on;
            xlabel('֡��'); ylabel('���');
            legend('λ����� (m)', '�Ƕ���� (��)');
            title('�������');
            
            subplot(obj.subplot_handles.stats);
            cla;
            mean_pos_error = mean(pos_error);
            std_pos_error = std(pos_error);
            mean_ang_error = mean(ang_error)*180/pi;
            std_ang_error = std(ang_error)*180/pi;
            text(0.1, 0.8, sprintf('λ�ù���:'), 'FontSize', 12);
            text(0.2, 0.7, sprintf('ƽ�����: %.2f m', mean_pos_error));
            text(0.2, 0.6, sprintf('��׼��: %.2f m', std_pos_error));
            text(0.1, 0.4, sprintf('�Ƕȹ���:'), 'FontSize', 12);
            text(0.2, 0.3, sprintf('ƽ�����: %.2f��', mean_ang_error));
            text(0.2, 0.2, sprintf('��׼��: %.2f��', std_ang_error));
            axis off;
            title('����ͳ��');
        end
        
        function show_final_results(obj, true_pos, est_pos, true_ang, est_ang, t_sim)
            figure('Name', '������������', 'Position', [150 150 1000 800]);
            subplot(2,2,1);
            pos_error = sqrt(sum((true_pos - est_pos).^2,2));
            ang_error = sqrt(sum((true_ang - est_ang).^2,2))*180/pi;
            [f_pos, x_pos] = ecdf(pos_error);
            [f_ang, x_ang] = ecdf(ang_error);
            plot(x_pos, f_pos, 'b-', 'LineWidth', 2);
            hold on;
            plot(x_ang, f_ang, 'r--', 'LineWidth', 2);
            grid on;
            xlabel('���'); ylabel('�ۻ�����');
            legend('λ����� (m)', '�Ƕ���� (��)');
            title('����ۻ��ֲ�');
            
            subplot(2,2,2);
            plot(t_sim, pos_error, 'b-', 'LineWidth', 1.5);
            hold on;
            plot(t_sim, ang_error, 'r--', 'LineWidth', 1.5);
            grid on;
            xlabel('ʱ�� (s)'); ylabel('���');
            legend('λ����� (m)', '�Ƕ���� (��)');
            title('���ʱ���ݻ�');
            
            subplot(2,2,[3,4]);
            rmse_pos = sqrt(mean(pos_error.^2));
            rmse_ang = sqrt(mean(ang_error.^2));
            max_pos_error = max(pos_error);
            max_ang_error = max(ang_error);
            text_info = {sprintf('λ�ù�������:'), sprintf('  RMSE: %.3f m', rmse_pos), sprintf('  ������: %.3f m', max_pos_error), sprintf('  ƽ�����: %.3f m', mean(pos_error)), sprintf('  ��׼��: %.3f m', std(pos_error)), '', sprintf('�Ƕȹ�������:'), sprintf('  RMSE: %.3f��', rmse_ang), sprintf('  ������: %.3f��', max_ang_error), sprintf('  ƽ�����: %.3f��', mean(ang_error)), sprintf('  ��׼��: %.3f��', std(ang_error))};
            text(0.1, 0.9, text_info, 'FontSize', 12, 'VerticalAlignment', 'top');
            axis off;
            title('��������ͳ��');
        end
    end
end

