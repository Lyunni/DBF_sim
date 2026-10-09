clc;clear;close all;

% ===== 物理常量 =====
c = 3e8;                % 光速 (m/s)

% ===== 仿真参数(可更改) =====
fc = 10e6;              % 载频 (Hz),可更改
lambda = c/fc;          % 工作波长 (m)

theta = 0;              % 波束指向俯仰角(相对+z轴),可更改,范围 0~pi
phi = 0;                % 波束指向方位角(相对+x轴),可更改,范围 0~2pi
theta_0 = 0;            % 输入信号到达方向俯仰角,可更改
phi_0 = 0;              % 输入信号到达方向方位角,可更改

% ===== 采样参数 =====
fs = 8 * fc;            % 采样率:每个载波周期采 8 个点
K = 128;                % K 个完整载波周期
N = K * (fs/fc);        % 采样点数(整周期数 -> FFT 无泄漏)
t = (0:N-1)/fs;         % 时间轴

% ===== 阵列参数(可更改) =====
array_m = 4;           % x 向阵元数,可更改
array_n = 12;            % y 向阵元数,可更改
array_dx = 12.8;        % x 向阵元间距 (m),可更改
array_dy = 12.8;        % y 向阵元间距 (m),可更改

% ===== 输入信号:参考阵元(0,0)收到的复载波 =====
omega = 2*pi*fc;        % 载波角频率
x_in = exp(1j*omega*t); % 复指数载波,幅度为 1

% ===== 波束形成:归一化阵列增益 =====
G = DBF_gain_2d(fc, theta, phi, array_m, array_n, array_dx, array_dy, theta_0, phi_0, x_in);
fprintf('归一化阵列增益 G = %.6f\n', G);

% ===== 二维方向图:固定波束指向(theta,phi),扫描信号来向(theta_0,phi_0) =====
theta_scan_deg = 0:1:90;        % 俯仰角扫描范围(°),可更改
phi_scan_deg = 0:1:360;         % 方位角扫描范围(°),可更改
theta_scan = deg2rad(theta_scan_deg);
phi_scan = deg2rad(phi_scan_deg);
[THETA_SCAN, PHI_SCAN] = meshgrid(theta_scan, phi_scan);

% 逐点扫描信号来向,每个角度调一次 DBF_gain_2d
G_pattern = zeros(size(THETA_SCAN));
for ii = 1:numel(theta_scan)
    for jj = 1:numel(phi_scan)
        G_pattern(jj, ii) = DBF_gain_2d(fc, theta, phi, array_m, array_n, ...
            array_dx, array_dy, THETA_SCAN(jj,ii), PHI_SCAN(jj,ii), x_in);
    end
end

floor_dB = -60;
G_pattern_dB = 20*log10(max(G_pattern, 10^(floor_dB/20)));

% ===== 图1: 三维波束方向图(球坐标) =====
% 用大窗口展示,真实数据比例,方向图实体占满坐标系
h2 = figure('Name', '3D Beam Pattern', 'NumberTitle', 'off', 'Position', [100 100 1200 1000]);

% 球坐标转换:半径=真实线性增益,theta=极角(相对+z),phi=方位角(相对+x)
% 不对数据做任何抬高或压低,直接使用 G_pattern
R = G_pattern;
X = R .* sin(THETA_SCAN) .* cos(PHI_SCAN);
Y = R .* sin(THETA_SCAN) .* sin(PHI_SCAN);
Z = R .* cos(THETA_SCAN);

% 画方向图表面,用 dB 值着色
s = surf(X, Y, Z, G_pattern_dB);
s.EdgeColor = [0.3 0.3 0.3];    % 深灰色网格边
s.FaceAlpha = 0.9;              % 轻微透明
s.LineWidth = 0.2;
hold on;

% 真实数据范围(不添加额外边距,让方向图占满坐标系)
xmax = max(abs(X(:)));
ymax = max(abs(Y(:)));
zmax = max(Z(:));
zmin = min(Z(:));

% 以 z 方向最大范围为基准设置立方体坐标系,让方向图实体占满画面 50% 以上
% 真实数据不变,只调整显示范围
box_limit = max([xmax, ymax, zmax, abs(zmin)]);
axis([-box_limit box_limit -box_limit box_limit zmin box_limit] * 1.05);
axis equal;

% 光照和材质
lighting gouraud;
camlight('headlight');
camlight('right');
material('dull');

grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('三维波束方向图 %d×%d 阵元, 指向 (θ=%.1f°, φ=%.1f°)', array_m, array_n, theta*180/pi, phi*180/pi));
view(120, 25);                  % 方位角 120°,仰角 25°

% 颜色条显示 dB
colorbar;
clim([floor_dB 0]);
colormap jet;
print(h2, '-dpng', '-r200', fullfile(fileparts(mfilename('fullpath')), 'DBF_2d_pattern_3d.png'));

% ===== 图2: 主切面方向图(phi=0°/90°/180°/270°) =====
h3 = figure('Name', 'Beam Pattern Cuts', 'NumberTitle', 'off');
[~, idx_phi0]   = min(abs(phi_scan_deg - 0));
[~, idx_phi90]  = min(abs(phi_scan_deg - 90));
[~, idx_phi180] = min(abs(phi_scan_deg - 180));
[~, idx_phi270] = min(abs(phi_scan_deg - 270));
plot(theta_scan_deg, G_pattern_dB(idx_phi0, :),   'LineWidth', 1.5); hold on;
plot(theta_scan_deg, G_pattern_dB(idx_phi90, :),  'LineWidth', 1.5);
plot(theta_scan_deg, G_pattern_dB(idx_phi180, :), 'LineWidth', 1.5);
plot(theta_scan_deg, G_pattern_dB(idx_phi270, :), 'LineWidth', 1.5);
xlabel('入射俯仰角 \theta_0 (°)');
ylabel('归一化增益 (dB)');
legend('\phi_0=0°', '\phi_0=90°', '\phi_0=180°', '\phi_0=270°', 'Location', 'best');
title('方向图主切面');
grid on;
ylim([floor_dB 0.5]);
print(h3, '-dpng', '-r200', fullfile(fileparts(mfilename('fullpath')), 'DBF_2d_pattern_cuts.png'));

% ===== 图3: 方位主平面方向图 =====
% theta 接近 0(近天顶指向)时,等 theta 圆锥整体落在主瓣 -3dB 之内,
% 方位切面退化(增益几乎不随 phi_0 变化,hpbw 求不到 -3dB 点),
% 此时改用 phi+90° 平面内过极点的正交切面代替。
hpbw_est_deg = rad2deg(0.886*lambda/min(array_m*array_dx, array_n*array_dy)); % 3dB 宽度粗估(°)
az_degenerate = rad2deg(theta) < hpbw_est_deg/2;

if ~az_degenerate
    % 常规:theta_0 = theta 固定,扫描 phi_0(0.1° 步进保证 3dB 宽度插值精度)
    az_axis_deg = 0:0.1:360;
    theta_az = theta * ones(size(az_axis_deg));
    phi_az = deg2rad(az_axis_deg);
    az_xlabel = '入射方位角 \phi_0 (°)';
    az_title = sprintf('方位主平面 (\\theta_0 = %.1f°)', rad2deg(theta));
else
    % 退化:phi+90° 平面内过峰带符号扫描(负角越过极点,等价于 |角| + 方位翻转 180°)
    az_ang_deg = rad2deg(theta) + (-15:0.1:15);
    az_axis_deg = az_ang_deg;
    theta_az = deg2rad(min(abs(az_ang_deg), 90));
    phi_az = phi + pi/2 + pi*(az_ang_deg < 0);
    az_xlabel = '相对波束指向的偏角 (°)';
    az_title = sprintf('正交切面 (近天顶指向, \\phi_0 = %.1f° 平面)', rad2deg(phi) + 90);
end

G_cut_az = zeros(size(az_axis_deg));     % 方位主平面增益
for jj = 1:numel(az_axis_deg)
    G_cut_az(jj) = DBF_gain_2d(fc, theta, phi, array_m, array_n, ...
        array_dx, array_dy, theta_az(jj), phi_az(jj), x_in);
end
G_cut_az_dB = 20*log10(max(G_cut_az, 10^(floor_dB/20)));

h4 = figure('Name', 'Azimuth Principal Plane Pattern', 'NumberTitle', 'off');
plot(az_axis_deg, G_cut_az_dB, 'LineWidth', 1.5); hold on;
if ~az_degenerate
    xline(rad2deg(phi), '--r', sprintf('波束指向 \\phi=%.1f°', rad2deg(phi)));
else
    xline(0, '--r', '波束指向');
end
xlabel(az_xlabel);
ylabel('归一化增益 (dB)');
title(az_title);
grid on; ylim([floor_dB 0.5]); xlim([az_axis_deg(1) az_axis_deg(end)]);
print(h4, '-dpng', '-r200', fullfile(fileparts(mfilename('fullpath')), 'DBF_2d_pattern_az_plane.png'));

% ===== 主瓣 3dB 宽度:过波束指向的全部取向切面,输出最大值/最小值 =====
% 三维方向图的主瓣截面一般不是圆(非方阵列各向异性更明显),3dB 宽度随切面
% 取向变化,单一主平面不能代表。这里以波束指向为轴,在 0~180° 内旋转切面
% 取向角 psi(psi=0 为俯仰向,psi=90° 为方位向),沿每个取向的大圆做以波束
% 为中心的带符号偏角扫描;切点用方向余弦表示,天顶/端射/下半球自动覆盖,
% 无退化分支。逐取向用 hpbw 求 3dB 宽度(见 hpbw.m)后取最大值/最小值。
% 理论值沿完全相同的切点用闭式阵因子 |sin(M·u/2)/(M·sin(u/2))| 的 x/y 向
% 乘积计算,u→0 处 0/0 取极限 1(与 1D 版 abs(u)<1e-9 处理一致)。
% 前提:单主瓣(无栅瓣);扫描窗宽按 1/cos(theta) 自适应加宽以覆盖端射展宽。
% 注意:终端输出不做 TeX 解释,希腊字母直接写 Unicode 字符。

% --- 波束指向处的局部正交坐标系(任意指向均有定义) ---
b  = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)]; % 波束方向单位矢量
e1 = [cos(theta)*cos(phi); cos(theta)*sin(phi); -sin(theta)]; % 俯仰向切矢量
e2 = [-sin(phi); cos(phi); 0];                                % 方位向切矢量

% --- 带符号偏角扫描窗:端射时主瓣按 1/cos(theta) 展宽,窗口自适应加宽 ---
half_win_deg = max(15, min(45, 2*hpbw_est_deg/max(abs(cos(theta)), 0.2)));
alpha_deg = -half_win_deg : 0.1 : half_win_deg;
alpha = deg2rad(alpha_deg);

psi_scan_deg = 0:5:175;             % 切面取向角(psi 与 psi+180° 为同一切面)
bw_sim = zeros(size(psi_scan_deg));
bw_theory = zeros(size(psi_scan_deg));

for kk = 1:numel(psi_scan_deg)
    psi = deg2rad(psi_scan_deg(kk));
    e_psi = cos(psi)*e1 + sin(psi)*e2;              % 该取向的切面方向
    D = b*cos(alpha) + e_psi*sin(alpha);            % 3×N 切点方向矢量(单位球面)
    theta_cut = acos(min(max(D(3,:), -1), 1));      % 数值保护后转球坐标
    phi_cut = atan2(D(2,:), D(1,:));

    % 仿真:沿切点逐点调 DBF_gain_2d
    G_cut = zeros(size(alpha));
    for ii = 1:numel(alpha)
        G_cut(ii) = DBF_gain_2d(fc, theta, phi, array_m, array_n, ...
            array_dx, array_dy, theta_cut(ii), phi_cut(ii), x_in);
    end
    G_cut_dB = 20*log10(max(G_cut, 10^(floor_dB/20)));

    % 理论:闭式阵因子,方向余弦差直接给出残差相位斜率
    ux = 2*pi*array_dx*(D(1,:) - b(1))/lambda;
    uy = 2*pi*array_dy*(D(2,:) - b(2))/lambda;
    AFx = abs(sin(array_m*ux/2) ./ (array_m*sin(ux/2)));
    AFy = abs(sin(array_n*uy/2) ./ (array_n*sin(uy/2)));
    AFx(abs(ux) < 1e-9) = 1;
    AFy(abs(uy) < 1e-9) = 1;
    AF_dB = 20*log10(max(AFx .* AFy, 10^(floor_dB/20)));

    % 峰值固定在偏角 0(波束指向)处,hpbw 直接求该取向的 3dB 宽度
    bw_sim(kk) = hpbw(G_cut_dB, alpha_deg);
    bw_theory(kk) = hpbw(AF_dB, alpha_deg);
end

n_nan = sum(isnan(bw_sim));
if n_nan > 0
    warning('%d 个取向的主瓣 -3dB 点触及扫描边界(±%.1f°),该取向按缺失处理。', ...
        n_nan, half_win_deg);
end

% --- 控制台输出:最大值/最小值(max/min 自动忽略 NaN) ---
[bw_sim_max, k_max] = max(bw_sim);
[bw_sim_min, k_min] = min(bw_sim);
bw_th_max = max(bw_theory);
bw_th_min = min(bw_theory);
fprintf('\n===== 主瓣 3dB 宽度(过波束指向的全部取向切面) =====\n');
fprintf('仿真结果: 最大值/最小值 = %.2f°/%.2f° (取向 ψ = %.0f°/%.0f°)\n', ...
    bw_sim_max, bw_sim_min, psi_scan_deg(k_max), psi_scan_deg(k_min));
fprintf('理论值:   最大值/最小值 = %.2f°/%.2f°\n', bw_th_max, bw_th_min);
fprintf('偏差:     %+.2f°/%+.2f° (%.2f%%/%.2f%%)\n', ...
    bw_sim_max - bw_th_max, bw_sim_min - bw_th_min, ...
    abs(bw_sim_max - bw_th_max)/abs(bw_th_max)*100, ...
    abs(bw_sim_min - bw_th_min)/abs(bw_th_min)*100);

