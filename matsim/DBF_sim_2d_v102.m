clc;clear;close all;

% ===== 物理常量 =====
c = 3e8;                % 光速 (m/s)

% ===== 仿真参数(可更改) =====
fc = 10e6;              % 载频 (Hz),可更改
lambda = c/fc;          % 工作波长 (m)

theta = pi/4;              % 波束指向俯仰角(相对+z轴),可更改,范围 0~pi
phi = pi/2;                % 波束指向方位角(相对+x轴),可更改,范围 0~2pi
theta_0 = 0;            % 输入信号到达方向俯仰角,可更改
phi_0 = 0;              % 输入信号到达方向方位角,可更改

% ===== 采样参数 =====
fs = 8 * fc;            % 采样率:每个载波周期采 8 个点
K = 128;                % K 个完整载波周期
N = K * (fs/fc);        % 采样点数(整周期数 -> FFT 无泄漏)
t = (0:N-1)/fs;         % 时间轴

% ===== 阵列参数(可更改) =====
array_m = 8;            % x 向阵元数,可更改
array_n = 8;           % y 向阵元数,可更改
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

% ===== 主瓣 3dB 宽度:从仿真结果求算,与理论值对比 =====
% 口径与 1D 版一致:找峰值 -> 两侧找 -3dB 点 -> interp1 插值 -> 相减得宽度(见 hpbw.m)。
% 理论值取均匀阵列阵因子闭式解 |sin(M·u/2)/(M·sin(u/2))| 的 x/y 向乘积,
% 沿与仿真完全相同的切点计算(退化分支同样是过极点带符号切面),
% u→0 处 0/0 取极限 1(与 1D 版 abs(u)<1e-9 处理一致)。
% 前提:单主瓣(无栅瓣)、非端射;越界情形不做一般化处理。
% 注意:终端输出不做 TeX 解释,希腊字母直接写 Unicode 字符。

% --- 精细俯仰切面(phi_0 = phi 所在平面) ---
% 近天顶指向时主瓣跨过 theta_0 = 0 边界,改用带符号扫描(负角越过极点,
% 等价于 |角| + 方位翻转 180°),与图3 退化分支同一处理。
if ~az_degenerate
    el_axis_deg = max(0, rad2deg(theta)-15) : 0.1 : min(90, rad2deg(theta)+15);
    theta_el = deg2rad(el_axis_deg);
    phi_el = phi * ones(size(el_axis_deg));
    el_label = sprintf('俯仰主平面 (φ0 = %.1f° 固定, 扫描 θ0)', rad2deg(phi));
else
    el_axis_deg = -15:0.1:15;                       % 相对波束指向的带符号偏角
    el_ang_deg = rad2deg(theta) + el_axis_deg;
    theta_el = deg2rad(abs(el_ang_deg));
    phi_el = phi + pi*(el_ang_deg < 0);
    el_label = sprintf('俯仰主平面 (近天顶指向, φ0 = %.1f° 平面过极点)', rad2deg(phi));
end

G_cut_el = zeros(size(el_axis_deg));
for ii = 1:numel(el_axis_deg)
    G_cut_el(ii) = DBF_gain_2d(fc, theta, phi, array_m, array_n, ...
        array_dx, array_dy, theta_el(ii), phi_el(ii), x_in);
end
G_cut_el_dB = 20*log10(max(G_cut_el, 10^(floor_dB/20)));

% --- 理论方向图(闭式阵因子,沿与仿真相同的切点) ---
% 残差相位斜率:ux = 2*pi*dx*(sin(theta_0)*cos(phi_0) - sin(theta)*cos(phi))/lambda
%             uy = 2*pi*dy*(sin(theta_0)*sin(phi_0) - sin(theta)*sin(phi))/lambda
ux_az = 2*pi*array_dx*(sin(theta_az).*cos(phi_az) - sin(theta)*cos(phi))/lambda;
uy_az = 2*pi*array_dy*(sin(theta_az).*sin(phi_az) - sin(theta)*sin(phi))/lambda;
AFx_az = abs(sin(array_m*ux_az/2) ./ (array_m*sin(ux_az/2)));
AFy_az = abs(sin(array_n*uy_az/2) ./ (array_n*sin(uy_az/2)));
AFx_az(abs(ux_az) < 1e-9) = 1;
AFy_az(abs(uy_az) < 1e-9) = 1;
AF_az_dB = 20*log10(max(AFx_az .* AFy_az, 10^(floor_dB/20)));

ux_el = 2*pi*array_dx*(sin(theta_el).*cos(phi_el) - sin(theta)*cos(phi))/lambda;
uy_el = 2*pi*array_dy*(sin(theta_el).*sin(phi_el) - sin(theta)*sin(phi))/lambda;
AFx_el = abs(sin(array_m*ux_el/2) ./ (array_m*sin(ux_el/2)));
AFy_el = abs(sin(array_n*uy_el/2) ./ (array_n*sin(uy_el/2)));
AFx_el(abs(ux_el) < 1e-9) = 1;
AFy_el(abs(uy_el) < 1e-9) = 1;
AF_el_dB = 20*log10(max(AFx_el .* AFy_el, 10^(floor_dB/20)));

% --- 求 3dB 宽度(相对各自切面峰值) ---
% 方位整圈扫描的峰值可能压在 0°/360° 边界上,求宽度前先循环平移使主瓣居中
% (宽度为两侧 -3dB 点之差,与平移量无关,不影响结果)
if ~az_degenerate
    [~, kp_az] = max(G_cut_az_dB);
    sh_az = round(numel(az_axis_deg)/2) - kp_az;
    G_az_c  = circshift(G_cut_az_dB, sh_az);
    AF_az_c = circshift(AF_az_dB,    sh_az);
else
    G_az_c  = G_cut_az_dB;
    AF_az_c = AF_az_dB;
end
bw_az_sim    = hpbw(G_az_c,  az_axis_deg);
bw_az_theory = hpbw(AF_az_c, az_axis_deg);
bw_el_sim    = hpbw(G_cut_el_dB, el_axis_deg);
bw_el_theory = hpbw(AF_el_dB,    el_axis_deg);

% --- 控制台输出 ---
if ~az_degenerate
    az_label = sprintf('方位主平面 (θ0 = %.1f° 固定, 扫描 φ0)', rad2deg(theta));
else
    az_label = sprintf('方位正交切面 (近天顶指向, φ0 = %.1f° 平面过极点)', rad2deg(phi)+90);
end
fprintf('\n===== 主瓣 3dB 宽度 =====\n');
fprintf('%s:\n', az_label);
fprintf('  仿真结果: %.2f°\n', bw_az_sim);
fprintf('  理论值:   %.2f°\n', bw_az_theory);
fprintf('  偏差:     %+.2f° (%.2f %%)\n', bw_az_sim - bw_az_theory, ...
    abs(bw_az_sim - bw_az_theory)/abs(bw_az_theory)*100);
fprintf('%s:\n', el_label);
fprintf('  仿真结果: %.2f°\n', bw_el_sim);
fprintf('  理论值:   %.2f°\n', bw_el_theory);
fprintf('  偏差:     %+.2f° (%.2f %%)\n', bw_el_sim - bw_el_theory, ...
    abs(bw_el_sim - bw_el_theory)/abs(bw_el_theory)*100);

