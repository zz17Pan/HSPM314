function pos = generate_true_trajectory(t)
    % ������ʵ�˶��켣
    % ����ʹ��һ���򵥵�Բ���˶�ʾ��
    %R = 40;     % �뾶
    %omega = 2*pi;  % ���ٶ�
    
    %x = R * cos(omega*t);
    x=6*t;
    %y = R * sin(omega*t);
    y=15*t;

    %z = 40 + 2*sin(2*omega*t);  % �߶�С������
    z=8*t;
    pos = [x, y, z];
end