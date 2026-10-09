clc


% ===== 红色对照(可更改) =====
fc0 = 12e6;             % 参考载频 (Hz):图2 红色虚线画 fc0 对应波长 λ0 的方向图
c = 3e8;
lamda0 = c/fc0;
% ===== 仿真参数(可更改) =====

fc = 10e6;              % 载频 (Hz),可更改
theta = 0;              % 天线主瓣方向(相对法线),可更改
array_num = 8;          % 阵元个数,可更改
theta_0 = pi/2;         % 输入信号到达方向(相对法线),可更改

% ===== 采样参数 =====
fs = 8 * fc;            % 采样率:每个载波周期采 8 个点
K = 128;                % K 个完整载波周期
N = K * (fs/fc);        % 采样点数(整周期数 -> FFT 无泄漏)
t = (0:N-1)/fs;         % 时间轴

% ===== 输入信号:参考阵元(0号)收到的复载波 =====
omega = 2*pi*fc;        % 载波角频率
phi = 0;                % 初相
x_in = exp(1j*(omega*t+phi)); % 复指数载波,幅度为 1

% ===== 波束形成:归一化阵列增益 =====
G = DBF_gain_1d(fc, theta, array_num, theta_0, x_in);
fprintf('归一化阵列增益 G = %.6f\n', G);

% ===== 天线方向图:固定主瓣指向 theta,扫描入射角 theta_0 =====
d = 12.8;                         % 阵元间距 (m),与 DBF_gain_1d 内部保持一致
lamda = 3e8/fc;                   % 工作波长 (m)
theta0_scan_deg = -90:0.1:90;     % 入射角扫描范围(相对法线,deg),0.1° 步进保证零点清晰
theta0_scan = deg2rad(theta0_scan_deg);

% 仿真链路(蓝色):主瓣指向固定为 theta,入射角逐点扫描,每个角度调一次 DBF_gain_1d
G_pattern = zeros(size(theta0_scan));
for k = 1:numel(theta0_scan)
    G_pattern(k) = DBF_gain_1d(fc, theta, array_num, theta0_scan(k), x_in);
end

% 理论方向图(闭式解,Dirichlet 核),u 为补偿后的残余相位步进;
% 图1 用当前波长 λ(与仿真逐点 golden 对比),图2 用参考波长 λ0(波长对照)
u0 = 2*pi*d*(sin(theta0_scan) - sin(theta))/lamda0;
AF_lamda0 = abs(sin(array_num*u0/2) ./ (array_num*sin(u0/2)));
AF_lamda0(abs(u0) < 1e-9) = 1;    % u→0 时 0/0 取极限 1(主瓣对准点)
u1 = 2*pi*d*(sin(theta0_scan) - sin(theta))/lamda;
AF_lamda = abs(sin(array_num*u1/2) ./ (array_num*sin(u1/2)));
AF_lamda(abs(u1) < 1e-9) = 1;

% 绘图用 dB(对数域:线性域看不出旁瓣;0 值压到 -60 dB 地板)
floor_dB = -60;
G_dB       = 20*log10(max(G_pattern,  10^(floor_dB/20)));
AF_lamda0_dB = 20*log10(max(AF_lamda0, 10^(floor_dB/20)));
AF_lamda_dB  = 20*log10(max(AF_lamda,  10^(floor_dB/20)));

% ===== 图1(第一张):同一波长 λ 下 理论 vs 仿真 对比(golden 验证)=====
figure;
plot(theta0_scan_deg, G_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta0_scan_deg, AF_lamda_dB, 'r--', 'LineWidth', 1.2);
grid on;
ylim([floor_dB, 0.5]);
xlabel('入射角 \theta_0 (^\circ)');
ylabel('归一化增益 (dB)');
title(sprintf(['理论 vs 仿真  N=%d, fc=%.1f MHz, d/\\lambda=%.3f'], ...
      array_num, fc/1e6, d/lamda));
legend('仿真 (DBF\_gain\_1d)', '理论 (Dirichlet)', 'Location', 'best');

% golden 数值核对:线性域最大偏差(两条曲线同式,只应差浮点精度)
fprintf('理论-仿真最大偏差 = %.3e\n', max(abs(G_pattern - AF_lamda)));

% ===== 图2(第二张):波长影响对比 蓝=当前 λ,红=参考 λ0 =====
figure;
plot(theta0_scan_deg, G_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta0_scan_deg, AF_lamda0_dB, 'r--', 'LineWidth', 1.2);
xline(rad2deg(theta), 'k--');    % 主瓣指向参考线
grid on;
ylim([floor_dB, 0.5]);
xlabel('入射角 \theta_0 (^\circ)');
ylabel('归一化增益 (dB)');
title(sprintf(['DBF 接收方向图  N=%d'],array_num));
legend(sprintf('当前 fc=%.1f MHz (\\lambda=%.1f m)', fc/1e6, lamda), ...
       sprintf('参考 fc0=%.1f MHz (\\lambda_0=%.1f m)', fc0/1e6, lamda0), ...
       '主瓣指向', 'Location', 'best');

% ===== 主瓣 3dB 宽度:直求,打印到终端 =====
% 前提:单主瓣(无栅瓣)、非端射;越界情形后续再处理。
% 注意:终端输出不做 TeX 解释,希腊字母直接写 Unicode 字符或 ASCII,不能用 \lambda
Pats = [G_dB; AF_lamda0_dB];               % 行1=蓝(当前λ),行2=红(λ0),峰值均 0 dB
for i = 1:2
    bw3(i) = hpbw(Pats(i,:), theta0_scan_deg);   % 找峰->两侧-3dB->插值,见 hpbw.m
end
display(['主瓣 3dB 宽度: ', num2str(bw3(1), '%.2f'), '° (当前 λ= ', num2str(lamda, '%.1f'), ' m)']);
display(['参考 3dB 宽度: ', num2str(bw3(2), '%.2f'), '° (参考 λ0= ', num2str(lamda0, '%.1f'), ' m)']);
