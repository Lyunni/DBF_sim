clc;clear;close all;

% ===== 物理常量 =====
c = 3e8;                % 光速 (m/s)

% ===== 仿真参数(可更改) =====
fc = 10e6;              % 载频 (Hz),可更改
lambda = c/fc;          % 工作波长 (m)

theta = pi/4;              % 波束指向俯仰角(相对+z轴),可更改,范围 0~pi
phi = pi/4;                % 波束指向方位角(相对+x轴),可更改,范围 0~2pi
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

% ===== 图3: 方位主平面方向图(固定 theta_0 = 波束俯仰角,扫描方位角) =====
phi_cut_deg = 0:0.5:360;                 % 方位主平面扫描范围(°)
phi_cut = deg2rad(phi_cut_deg);

G_cut_az = zeros(size(phi_cut));         % 方位主平面增益
for jj = 1:numel(phi_cut)
    G_cut_az(jj) = DBF_gain_2d(fc, theta, phi, array_m, array_n, ...
        array_dx, array_dy, theta, phi_cut(jj), x_in);
end
G_cut_az_dB = 20*log10(max(G_cut_az, 10^(floor_dB/20)));

h4 = figure('Name', 'Azimuth Principal Plane Pattern', 'NumberTitle', 'off');
plot(phi_cut_deg, G_cut_az_dB, 'LineWidth', 1.5); hold on;
xline(rad2deg(phi), '--r', sprintf('波束指向 \\phi=%.1f°', rad2deg(phi)));
xlabel('入射方位角 \phi_0 (°)');
ylabel('归一化增益 (dB)');
title(sprintf('方位主平面 (\\theta_0 = %.1f°)', rad2deg(theta)));
grid on; ylim([floor_dB 0.5]); xlim([phi_cut_deg(1) phi_cut_deg(end)]);
print(h4, '-dpng', '-r200', fullfile(fileparts(mfilename('fullpath')), 'DBF_2d_pattern_az_plane.png'));
