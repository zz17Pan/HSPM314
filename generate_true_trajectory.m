function pos = generate_true_trajectory(t)
    % 生成真实运动轨迹
    % 这里使用一个简单的圆周运动示例
    %R = 40;     % 半径
    %omega = 2*pi;  % 角速度
    
    %x = R * cos(omega*t);
    x=6*t;
    %y = R * sin(omega*t);
    y=15*t;

    %z = 40 + 2*sin(2*omega*t);  % 高度小幅波动
    z=8*t;
    pos = [x, y, z];
end