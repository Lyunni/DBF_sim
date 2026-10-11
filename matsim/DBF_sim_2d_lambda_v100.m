clc;clear;close all;

% DBF_sim_2d_lambda_v100
% 以 DBF_sim_2d_v102 为蓝本,研究信号固定沿主瓣方向入射时,
% 阵列增益与信号波长的关系,并给出栅瓣条件的边界值。
%
% 要点:蓝本中 DBF 补偿权按信号载频 fc 实时计算(权随频率自适应),
% 信号沿主瓣方向入射时增益恒为 1,看不到频率依赖。实际相控阵的
% 移相器按设计波长 lambda_ref 装定后固定不变,信号波长偏离设计值时
% 主瓣偏移(波束斜视),信号方向增益下降;波长小到满足栅瓣条件时
% 可见空间内出现与主瓣等高的栅瓣。本脚本两种情形都画出来对比。

% ===== 物理常量 =====
c = 3e8;                % 光速 (m/s)

% ===== 仿真参数(可更改) =====
fc_ref = 10e6;              % 设计载频 (Hz),补偿权按此频率装定
lambda_ref = c/fc_ref;      % 设计波长 (m)

theta = deg2rad(45);        % 主瓣(波束)指向俯仰角,可更改,范围 0~pi
                            % 注意:theta=0 时补偿权为全 1,增益恒为 1,无频率依赖
phi = 0;                    % 主瓣指向方位角,可更改,范围 0~2pi

% 信号固定沿主瓣方向入射
theta_0 = theta;
phi_0 = phi;

% ===== 阵列参数(与蓝本一致,可更改) =====
array_m = 8;            % x 向阵元数
array_n = 8;            % y 向阵元数
array_dx = 12.8;        % x 向阵元间距 (m)
array_dy = 12.8;        % y 向阵元间距 (m)

% ===== 波长扫描范围(可更改) =====
lambda_scan = 8:0.05:60;    % 信号波长扫描 (m),对应 fc = c/lambda 约 5~37.5 MHz

% ===== 方向余弦 =====
u0 = sin(theta)*cos(phi);   % 主瓣/信号方向 x 向方向余弦
v0 = sin(theta)*sin(phi);   % 主瓣/信号方向 y 向方向余弦

% ===== 增益与波长的关系(闭式阵因子) =====
% 信号沿主瓣方向入射,阵元(m,n)残差相位:
%   psi_x = 2*pi*dx*u0*(1/lambda - 1/lambda_ref)   (固定权,权按 lambda_ref 装定)
% 增益 = |AFx| .* |AFy|,与 DBF_gain_2d 的同相叠加幅度/(M*N) 严格相等
% (相干叠加的时间域最大值 = 复数和的模),u→0 处 0/0 取极限 1(同蓝本处理)。
psi_x = 2*pi*array_dx*u0*(1./lambda_scan - 1/lambda_ref);
psi_y = 2*pi*array_dy*v0*(1./lambda_scan - 1/lambda_ref);

AFx = abs(sin(array_m*psi_x/2) ./ (array_m*sin(psi_x/2)));
AFy = abs(sin(array_n*psi_y/2) ./ (array_n*sin(psi_y/2)));
AFx(abs(psi_x) < 1e-9) = 1;
AFy(abs(psi_y) < 1e-9) = 1;

G_fix = AFx .* AFy;                 % 固定权(按 lambda_ref 装定):增益随波长变化
G_adapt = ones(size(lambda_scan));  % 自适应权(蓝本,权随 fc 重算):恒为 1

floor_dB = -60;
G_fix_dB = 20*log10(max(G_fix, 10^(floor_dB/20)));

% ===== 用蓝本函数 DBF_gain_2d 抽查验证:权随频率自适应时增益恒为 1 =====
fprintf('===== 蓝本算法抽查(自适应权,信号沿主瓣入射) =====\n');
lambda_check = [12 20 30 45 60];
K = 32;                             % 载波周期数(抽查用,比蓝本 128 小以提速)
for lam = lambda_check
    fc = c/lam;
    fs = 8*fc;
    t = (0:K*fs/fc-1)/fs;
    x_in = exp(1j*2*pi*fc*t);
    G = DBF_gain_2d(fc, theta, phi, array_m, array_n, ...
        array_dx, array_dy, theta_0, phi_0, x_in);
    fprintf('  lambda = %5.1f m (fc = %6.2f MHz): G = %.6f\n', lam, fc/1e6, G);
end

% ===== 栅瓣条件边界值 =====
% 矩形阵列沿 x/y 两轴独立:方向余弦空间内栅瓣间距为 lambda/d,
% 主瓣位于 (u0,v0),最近的栅瓣进入可见区(|u|<=1)的边界为
%   lambda_bx = dx*(1+|u0|),  lambda_by = dy*(1+|v0|)
% lambda < lambda_b 时该轴向出现栅瓣;无栅瓣要求 lambda >= max(lambda_bx, lambda_by)。
% 特例:端射(theta=90°)边界 lambda = 2d,即全空域无栅瓣的经典条件 d <= lambda/2。
lam_bx = array_dx*(1 + abs(u0));
lam_by = array_dy*(1 + abs(v0));
lam_b = max(lam_bx, lam_by);

fprintf('\n===== 栅瓣条件边界值(指向 theta=%.1f°, phi=%.1f°) =====\n', ...
    rad2deg(theta), rad2deg(phi));
fprintf('方向余弦: u0 = %+.4f, v0 = %+.4f\n', u0, v0);
fprintf('x 向栅瓣边界: lambda_bx = dx*(1+|u0|) = %.2f m (fc > %.2f MHz 时出现 x 向栅瓣)\n', ...
    lam_bx, c/lam_bx/1e6);
fprintf('y 向栅瓣边界: lambda_by = dy*(1+|v0|) = %.2f m (fc > %.2f MHz 时出现 y 向栅瓣)\n', ...
    lam_by, c/lam_by/1e6);
fprintf('无栅瓣条件:   lambda >= %.2f m  (fc <= %.2f MHz)\n', lam_b, c/lam_b/1e6);
fprintf('全空域扫描无栅瓣(端射特例): lambda >= 2d = %.2f m, 即 d <= lambda/2\n', ...
    2*max(array_dx, array_dy));

% ===== 图1: 增益与信号波长的关系 =====
h1 = figure('Name', 'Gain vs Wavelength', 'NumberTitle', 'off', 'Position', [80 80 1100 600]);
plot(lambda_scan, G_fix_dB, 'LineWidth', 1.5); hold on;
plot(lambda_scan, 20*log10(G_adapt), '--', 'LineWidth', 1.5);
xline(lambda_ref, '-.', sprintf('设计波长 \\lambda_{ref} = %.1f m', lambda_ref), ...
    'LabelHorizontalAlignment', 'left');
xline(lam_bx, '--', sprintf('x 向栅瓣边界 %.2f m', lam_bx), ...
    'LabelHorizontalAlignment', 'left');
xline(lam_by, '--', sprintf('y 向栅瓣边界 %.2f m', lam_by), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
xlabel('信号波长 \lambda (m)');
ylabel('归一化增益 (dB)');
title(sprintf('增益-波长关系: %d×%d 阵元, 主瓣 (\\theta=%.1f°, \\phi=%.1f°), 信号沿主瓣入射', ...
    array_m, array_n, rad2deg(theta), rad2deg(phi)));
legend('固定权(按 \lambda_{ref} 装定)', '自适应权(蓝本, 权随频率重算)', 'Location', 'best');
grid on;
ylim([floor_dB 3]);
xlim([lambda_scan(1) lambda_scan(end)]);
print(h1, '-dpng', '-r200', fullfile(fileparts(mfilename('fullpath')), 'DBF_2d_lambda_gain.png'));

% ===== 图2: 栅瓣区内的方向图切面(固定权) =====
% 取一个小于栅瓣边界的演示波长,在 phi 切面内做带符号偏角扫描,
% 可见主瓣偏离信号方向(波束斜视)且出现与主瓣等高的栅瓣。
lambda_demo = max(lambda_scan(1), min(lambda_scan(end), round(0.7*lam_b, 1)));
alpha_deg = -90:0.1:90;
alpha = deg2rad(alpha_deg);
u_cut = sin(alpha)*cos(phi);        % phi 切面内的 x 向方向余弦
v_cut = sin(alpha)*sin(phi);        % phi 切面内的 y 向方向余弦

% 固定权:残差 = 2*pi*d*(u/lambda_demo - u0/lambda_ref)
px_fix = 2*pi*array_dx*(u_cut/lambda_demo - u0/lambda_ref);
py_fix = 2*pi*array_dy*(v_cut/lambda_demo - v0/lambda_ref);
% 自适应权(蓝本):残差 = 2*pi*d*(u - u0)/lambda_demo
px_ad = 2*pi*array_dx*(u_cut - u0)/lambda_demo;
py_ad = 2*pi*array_dy*(v_cut - v0)/lambda_demo;

AFx_fix = abs(sin(array_m*px_fix/2) ./ (array_m*sin(px_fix/2)));
AFy_fix = abs(sin(array_n*py_fix/2) ./ (array_n*sin(py_fix/2)));
AFx_fix(abs(px_fix) < 1e-9) = 1;
AFy_fix(abs(py_fix) < 1e-9) = 1;
AFx_ad = abs(sin(array_m*px_ad/2) ./ (array_m*sin(px_ad/2)));
AFy_ad = abs(sin(array_n*py_ad/2) ./ (array_n*sin(py_ad/2)));
AFx_ad(abs(px_ad) < 1e-9) = 1;
AFy_ad(abs(py_ad) < 1e-9) = 1;

G_cut_fix_dB = 20*log10(max(AFx_fix .* AFy_fix, 10^(floor_dB/20)));
G_cut_ad_dB  = 20*log10(max(AFx_ad  .* AFy_ad,  10^(floor_dB/20)));

% 固定权主瓣位置:sin(alpha_ml) = sin(theta)*lambda_demo/lambda_ref
sin_ml = sin(theta)*lambda_demo/lambda_ref;

h2 = figure('Name', 'Grating Lobe Pattern Cut', 'NumberTitle', 'off', 'Position', [80 80 1100 600]);
plot(alpha_deg, G_cut_fix_dB, 'LineWidth', 1.5); hold on;
plot(alpha_deg, G_cut_ad_dB, '--', 'LineWidth', 1.5);
xline(rad2deg(theta), '-.', sprintf('信号/主瓣指向 %.1f°', rad2deg(theta)));
if abs(sin_ml) <= 1
    xline(rad2deg(asin(sin_ml)), ':', sprintf('固定权主瓣偏移至 %.1f°', rad2deg(asin(sin_ml))));
end
xlabel(sprintf('\\phi = %.1f° 切面内的带符号俯仰偏角 (°)', rad2deg(phi)));
ylabel('归一化增益 (dB)');
title(sprintf('\\lambda = %.1f m < 栅瓣边界 %.2f m: 固定权方向图出现栅瓣', lambda_demo, lam_b));
legend('固定权(按 \lambda_{ref} 装定)', '自适应权(蓝本)', 'Location', 'best');
grid on;
ylim([floor_dB 3]);
xlim([alpha_deg(1) alpha_deg(end)]);
print(h2, '-dpng', '-r200', fullfile(fileparts(mfilename('fullpath')), 'DBF_2d_lambda_grating_cut.png'));

fprintf('\n演示波长 lambda_demo = %.1f m(< 栅瓣边界 %.2f m),方向图切面见图2。\n', ...
    lambda_demo, lam_b);
