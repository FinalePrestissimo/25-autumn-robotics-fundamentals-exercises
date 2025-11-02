clear; clc; close all;

%% 参数设置
d1 = 0.3;  % m
d3 = 0.5;  % m

% 初始关节角
theta0_deg = [20; -50; 100];
theta0 = deg2rad(theta0_deg);

% 初始和终止姿态角（xyz动轴欧拉角）
phi0_deg = [-157.8240; 46.0418; 70.4798];
phif_deg = [174.9616; 8.6492; 90.3813];
phi0 = deg2rad(phi0_deg);
phif = deg2rad(phif_deg);

% 时间参数
tf = 10;         % 总时间 10s
dt = 0.1;        % 采样周期 0.1s
t = 0:dt:tf;
N = length(t);



%% 方法1
tau1 = t / tf;
s1 = 3*tau1.^2 - 2*tau1.^3;
s1_dot = (6*tau1 - 6*tau1.^2) / tf;

phi1_traj = phi0 + (phif - phi0) .* s1;
phi1_dot_traj = (phif - phi0) .* s1_dot;
phi1_traj_deg = rad2deg(phi1_traj);



%% 方法2）
R0 = euler_xyz_to_rotation(phi0);
Rf = euler_xyz_to_rotation(phif);
R_rel = R0' * Rf;
[theta_eq, k_axis] = rotation_to_axis_angle(R_rel);

tau2 = t / tf;
s2 = 10*tau2.^3 - 15*tau2.^4 + 6*tau2.^5;
theta_traj = theta_eq * s2;

phi2_traj = zeros(3, N);
phi2_dot_traj = zeros(3, N);

for i = 1:N
    R_current = axis_angle_to_rotation(k_axis, theta_traj(i));
    R_abs = R0 * R_current;
    phi2_traj(:, i) = rotation_to_euler_xyz(R_abs);
    
    if i > 1
        phi2_dot_traj(:, i) = (phi2_traj(:, i) - phi2_traj(:, i-1)) / dt;
    end
end
phi2_dot_traj(:, 1) = phi2_dot_traj(:, 2);
phi2_traj_deg = rad2deg(phi2_traj);



%% 计算关节角轨迹
theta1_traj = zeros(3, N);
for i = 1:N
    if i == 1
        theta_init = theta0;
    else
        theta_init = theta1_traj(:, i-1);
    end
    
    R_target = euler_xyz_to_rotation(phi1_traj(:, i));
    theta1_traj(:, i) = inverse_kinematics_wrist(R_target, d1, d3, theta_init);
end
theta1_traj_deg = rad2deg(theta1_traj);

theta2_traj = zeros(3, N);
for i = 1:N
    if i == 1
        theta_init = theta0;
    else
        theta_init = theta2_traj(:, i-1);
    end
    
    R_target = euler_xyz_to_rotation(phi2_traj(:, i));
    theta2_traj(:, i) = inverse_kinematics_wrist(R_target, d1, d3, theta_init);
end
theta2_traj_deg = rad2deg(theta2_traj);



%% 绘图
figure('Name', '末端姿态角变化曲线');
euler_labels = {'\phi (roll)', '\theta (pitch)', '\psi (yaw)'};
for i = 1:3
    subplot(3, 1, i);
    plot(t, phi1_traj_deg(i, :), 'b-', 'LineWidth', 1.5, 'DisplayName', '方法1');
    hold on;
    plot(t, phi2_traj_deg(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', '方法2');
    grid on;
    xlabel('时间 (s)');
    ylabel([euler_labels{i} ' (度)']);
    legend('Location', 'best');
end

figure('Name', '末端角速度曲线');
for i = 1:3
    subplot(3, 1, i);
    plot(t, phi1_dot_traj(i, :), 'b-', 'LineWidth', 1.5, 'DisplayName', '方法1');
    hold on;
    plot(t, phi2_dot_traj(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', '方法2');
    grid on;
    xlabel('时间 (s)');
    ylabel(['\omega_' num2str(i) ' (rad/s)']);
    legend('Location', 'best');
end

figure('Name', '关节角曲线');
for i = 1:3
    subplot(3, 1, i);
    plot(t, theta1_traj_deg(i, :), 'b-', 'LineWidth', 1.5, 'DisplayName', '方法1');
    hold on;
    plot(t, theta2_traj_deg(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', '方法2');
    grid on;
    xlabel('时间 (s)');
    ylabel(['\theta_' num2str(i) ' (度)']);
    legend('Location', 'best');
end

%% 辅助函数

% xyz动轴欧拉角转旋转矩阵（按 Rx * Ry * Rz 顺序）
function R = euler_xyz_to_rotation(phi)
    alpha = phi(1);  % 绕x轴
    beta = phi(2);   % 绕y轴
    gamma = phi(3);  % 绕z轴
    
    Rx = [1, 0, 0;
          0, cos(alpha), -sin(alpha);
          0, sin(alpha), cos(alpha)];
    
    Ry = [cos(beta), 0, sin(beta);
          0, 1, 0;
          -sin(beta), 0, cos(beta)];
    
    Rz = [cos(gamma), -sin(gamma), 0;
          sin(gamma), cos(gamma), 0;
          0, 0, 1];
    
    % 动轴欧拉角：先绕固定轴X，再绕新Y，最后绕新Z
    R = Rx * Ry * Rz;
end

% 旋转矩阵转xyz动轴欧拉角
function phi = rotation_to_euler_xyz(R)
    % 从 R = Rx * Ry * Rz 提取欧拉角
    % R(3,2) = sin(alpha)*cos(beta)*sin(gamma) + cos(alpha)*cos(gamma)
    % R(2,3) = sin(alpha)*sin(gamma) + cos(alpha)*cos(beta)*cos(gamma)
    
    beta = asin(R(1,3));  % sin(beta) = R(1,3)
    
    if abs(cos(beta)) > 1e-6
        alpha = atan2(-R(2,3), R(3,3));
        gamma = atan2(-R(1,2), R(1,1));
    else
        % 万向锁情况
        alpha = 0;
        gamma = atan2(R(2,1), R(2,2));
    end
    
    phi = [alpha; beta; gamma];
end

% 旋转矩阵转等效轴角
function [theta, k] = rotation_to_axis_angle(R)
    theta = acos((trace(R) - 1) / 2);
    
    if abs(theta) < 1e-6
        k = [0; 0; 1];
        theta = 0;
    else
        k = (1 / (2*sin(theta))) * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
        k = k / norm(k);
    end
end

% 轴角转旋转矩阵（罗德里格斯公式）
function R = axis_angle_to_rotation(k, theta)
    K = [0, -k(3), k(2);
         k(3), 0, -k(1);
         -k(2), k(1), 0];
    
    R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);
end

% 球腕逆运动学
function theta = inverse_kinematics_wrist(R_target, d1, d3, theta_init)
    % 简化的球腕逆运动学（基于数值优化）
    objective = @(th) norm(compute_wrist_rotation(th, d1, d3) - R_target, 'fro')^2;
    
    options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 500);
    theta = fminsearch(objective, theta_init, options);
end

% 计算球腕旋转矩阵
function R = compute_wrist_rotation(theta, d1, d3)
    % 3R球腕的旋转矩阵计算
    c1 = cos(theta(1)); s1 = sin(theta(1));
    c2 = cos(theta(2)); s2 = sin(theta(2));
    c3 = cos(theta(3)); s3 = sin(theta(3));
    
    % 简化的球腕旋转
    R = [c1*c2*c3-s1*s3, -c1*c2*s3-s1*c3, c1*s2;
         s1*c2*c3+c1*s3, -s1*c2*s3+c1*c3, s1*s2;
         -s2*c3, s2*s3, c2];
end
