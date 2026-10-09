function bw = hpbw(P, axis_deg)
%HPBW 从 dB 方向图曲线求主瓣 3dB 宽度(°)
%   P 为 dB 曲线(峰值 0 dB),axis_deg 为对应的角度轴(°)。
%   某侧找不到 -3dB 点(主瓣触及扫描边界)时返回 NaN 并打印警告。
[pk, kp] = max(P);                       % 主瓣中心
th = pk - 3;                             % -3dB 阈值(相对本曲线峰值)
il = find(P(1:kp)   <= th, 1, 'last');   % 左 -3dB 所在区间下沿
ir = kp - 1 + find(P(kp:end) <= th, 1);  % 右 -3dB 所在区间上沿
if isempty(il) || isempty(ir)
    warning('主瓣 -3dB 点触及扫描边界,无法求宽度,请扩大扫描范围。');
    bw = NaN;
    return;
end
th_L = interp1(P(il:il+1), axis_deg(il:il+1), th);
th_R = interp1(P(ir-1:ir), axis_deg(ir-1:ir), th);
bw = th_R - th_L;
end
