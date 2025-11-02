clear; clc; close all;

%% 参数设置
% 机械臂DH参数
d1 = 0.5;  % m
a2 = 0.4;  % m
a3 = 0.6;  % m

% 关节初值
theta0_deg = [26.5651; -126.9498; 87.6120];
theta0 = deg2rad(theta0_deg);  % 转换为弧度

% 关键点坐标
p0 = [0.2; 0.1; 1.2];      % 起点 P_0
pf = [0.1; -0.2; 0.8];     % 终点 P_f
Oc = [0; 0; 1];            % 圆心 O_c

% 时间参数
tf = 100;        % 总时间 100s
dt = 0.1;        % 采样周期 0.1s
t = 0:dt:tf;     % 时间序列
N = length(t);   % 采样点数



%% 圆弧轨迹规划
% 圆弧半径
r = norm(p0 - Oc);
fprintf('圆弧半径 r = %.4f m\n\n', r);

% 构建圆平面坐标系基向量 i, j, k
i_vec = (p0 - Oc) / norm(p0 - Oc);
fprintf('i = [%.4f, %.4f, %.4f]^T\n', i_vec(1), i_vec(2), i_vec(3));

v0 = p0 - Oc;
vf = pf - Oc;
n = cross(v0, vf);
fprintf('法向量 n = [%.4f, %.4f, %.4f]^T\n', n(1), n(2), n(3));

k_vec = n / norm(n);
fprintf('k = [%.4f, %.4f, %.4f]^T\n', k_vec(1), k_vec(2), k_vec(3));

j_vec = cross(k_vec, i_vec);
fprintf('j = [%.4f, %.4f, %.4f]^T\n', j_vec(1), j_vec(2), j_vec(3));

% 从 P_0 到 P_f 的圆心角
cos_phi_f = dot(v0, vf) / (norm(v0) * norm(vf));
phi_f = acos(cos_phi_f);
fprintf('P_0 到 P_f 的圆心角 phi_f = %.2f度\n', rad2deg(phi_f));



%% 时间规划（三次多项式）
tau = t / tf;  % 归一化时间 [0,1]
lambda = 3*tau.^2 - 2*tau.^3;
lambda_dot = (6*tau - 6*tau.^2) / tf;
lambda_ddot = (6 - 12*tau) / tf^2;



%% 生成圆弧轨迹
p_traj = zeros(3, N);

phi_0 = 0;  % 起点对应角度为 0
for i = 1:N
    % 当前角度
    phi = phi_0 + lambda(i) * (phi_f - phi_0);
    % 圆弧轨迹
    p_traj(:, i) = Oc + r * (cos(phi) * i_vec + sin(phi) * j_vec);
end



%% 逆运动学求解
theta_traj = zeros(3, N);

% 定义正运动学误差函数
fk_error = @(theta, p_target) forward_kinematics(theta, d1, a2, a3) - p_target;

% 设置求解器选项
options = optimoptions('fsolve', ...
    'Display', 'off', ...           % 不显示迭代信息
    'TolFun', 1e-8, ...             % 函数容差
    'TolX', 1e-8, ...               % 变量容差
    'MaxIterations', 1000);         % 最大迭代次数

fprintf('\n开始逆运动学求解...\n');

for i = 1:N
    p_target = p_traj(:, i);
    
    % 使用前一时刻的解作为初值（第一次使用给定初值）
    if i == 1
        theta_init = theta0;
    else
        theta_init = theta_traj(:, i-1);
    end
    
    % 使用fsolve求解
    [theta_sol, fval, exitflag] = fsolve(@(theta) fk_error(theta, p_target), ...
                                          theta_init, options);
    
    % 检查求解状态
    if exitflag <= 0
        warning('在时刻 t=%.2f 处逆运动学求解未收敛，exitflag=%d', t(i), exitflag);
    end
    
    theta_traj(:, i) = theta_sol;
    
    % 每100个点显示一次进度
    if mod(i, 100) == 0
        fprintf('已完成 %d/%d 点 (%.1f%%)\n', i, N, 100*i/N);
    end
end

fprintf('逆运动学求解完成！\n');

% 转换为角度
theta_traj_deg = rad2deg(theta_traj);



%% 验证正运动学
% 验证起点的正运动学
p0_verify = forward_kinematics(theta0, d1, a2, a3);
fprintf('\n正运动学验证（初始关节角对应位置）：\n');
fprintf('关节角: theta = [%.4f, %.4f, %.4f]^T (度)\n', theta0_deg(1), theta0_deg(2), theta0_deg(3));
fprintf('给定起点: P_0 = [%.4f, %.4f, %.4f]^T\n', p0(1), p0(2), p0(3));
fprintf('正运动学计算: P = [%.4f, %.4f, %.4f]^T\n', p0_verify(1), p0_verify(2), p0_verify(3));
fprintf('误差: %.6f m\n', norm(p0 - p0_verify));



%% 绘图
% 关节角曲线
figure('Name', '关节角曲线', 'Position', [100, 100, 1200, 800]);
for i = 1:3
    subplot(3, 1, i);
    plot(t, theta_traj_deg(i, :), 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('时间 (s)', 'FontSize', 12);
    ylabel(['\theta_' num2str(i) ' (度)'], 'FontSize', 12);
    title(['关节 ' num2str(i) ' 角度曲线'], 'FontSize', 14);
    xlim([0, tf]);
end

% 末端位置曲线
figure('Name', '末端位置曲线', 'Position', [150, 150, 1200, 800]);
coords = {'x', 'y', 'z'};
for i = 1:3
    subplot(3, 1, i);
    plot(t, p_traj(i, :), 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('时间 (s)', 'FontSize', 12);
    ylabel([coords{i} ' (m)'], 'FontSize', 12);
    title(['末端' coords{i} '坐标曲线'], 'FontSize', 14);
    xlim([0, tf]);
end

% 3D轨迹
figure('Name', '3D轨迹', 'Position', [200, 200, 800, 800]);
plot3(p_traj(1, :), p_traj(2, :), p_traj(3, :), 'b-', 'LineWidth', 2);
hold on;

% 标记关键点
plot3(p0(1), p0(2), p0(3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'P_0（起点）');
plot3(pf(1), pf(2), pf(3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'P_f（终点）');
plot3(Oc(1), Oc(2), Oc(3), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'O_c（圆心）');

% 绘制从圆心到关键点的连线
plot3([Oc(1), p0(1)], [Oc(2), p0(2)], [Oc(3), p0(3)], 'g--', 'LineWidth', 1);
plot3([Oc(1), pf(1)], [Oc(2), pf(2)], [Oc(3), pf(3)], 'r--', 'LineWidth', 1);

grid on;
xlabel('x (m)', 'FontSize', 12);
ylabel('y (m)', 'FontSize', 12);
zlabel('z (m)', 'FontSize', 12);
title('末端圆弧轨迹（3D视图）', 'FontSize', 14);
legend('Location', 'best');
axis equal;
view(45, 30);



%% 辅助函数
% 正运动学函数
function p = forward_kinematics(theta, d1, a2, a3)
    theta1 = theta(1);
    theta2 = theta(2);
    theta3 = theta(3);
    
    x = cos(theta1) * (a2*cos(theta2) + a3*cos(theta2+theta3));
    y = sin(theta1) * (a2*cos(theta2) + a3*cos(theta2+theta3));
    z = d1 - a2*sin(theta2) - a3*sin(theta2+theta3);
    
    p = [x; y; z];
end
