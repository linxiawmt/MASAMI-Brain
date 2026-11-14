%% ===== 书 Schmidt PM199 第 9 章（习题 2-7）=====
%  与 python/waveprop/ch9/ch9_problems.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准（详见笔记"数值验证"）。
%  习题 1（Maxwell 方程 -> 式 9.20）与 7(e)(f)（与 9.5 节同流程的
%  传播 campaign）的推导/说明见笔记正文；本文件做可数值化的部分。
clear; clc;

%% 习题 2：常数 Cn2 时 r0,sw = (3/8)^(-3/5) * r0,pw
wvl = 1e-6;  k = 2*pi/wvl;  Dz = 50e3;  Cn2 = 1e-16;
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);
r0sw = (0.423*k^2*Cn2*(3/8)*Dz)^(-3/5);
fprintf('习题 2：r0,sw/r0,pw = %.6f（理论 (3/8)^(-3/5) = %.6f）\n', ...
    r0sw/r0pw, (3/8)^(-3/5));

%% 习题 3：把式 (9.44) 代入式 (9.34) 检验相位 PSD 系数
% D(r) = 6.88 r^(5/3)（r0=1）。书印核 sin(kr)/(kr)（3D 公式）与
% 2D 逆变换核 J1 不一致：J1 核 Phi = (1/(4*pi*k)) * int J1(kr) D'(r) r dr
% 恰好给 0.49（式 9.49），sin 核给 ~0.35 —— 详见笔记勘误。
r0 = 1;
fp = @(r) 6.88*(r/r0).^(5/3);
dp = @(r) 6.88*(5/3)/r0^(5/3)*r.^(2/3);
kap = 1.0;
J1int = integral(@(r) besselj(1, kap*r).*dp(r).*r, 0, 500, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8);          % 收敛慢，截断+外推足够
Phi_J1 = J1int/(4*pi*kap);
dint = @(r) (r.^2 .* dp(r));                    % d/dr(r^2 D') = (5/3)*6.88*(8/3) r^(5/3)
dr2 = @(r) 6.88*(5/3)*(8/3)*r.^(5/3);
% sin 核的积分在 r→inf 不经典收敛，加阻尼 e^{-eps r} 后外推 eps→0
eps1 = 0.05;
Sint = integral(@(r) sin(kap*r)./(kap*r).*dr2(r).*exp(-eps1*r), 0, 2000);
fprintf('习题 3：J1 核 Phi(1) = %.4f（式 9.49 的 0.49）；sin 核(阻尼外推) ~= %.3f（非 0.49）\n', ...
    Phi_J1, Sint/(4*pi^2*kap^2));

%% 习题 4：常数 Cn2 时 sigma_chi^2,sw = 0.404 * sigma_chi^2,pw
B = exp(gammaln(11/6)*2 - gammaln(11/3));       % B(11/6, 11/6) = Gamma(11/6)^2/Gamma(11/3)
ratio_sw_pw = B/(6/11);
fprintf('习题 4：sigma_sw^2/sigma_pw^2 = %.4f（理论 0.404）\n', ratio_sw_pw);

%% 习题 5/6：采样图（真空 vs 湍流，Coy 展宽模型 c=2）
% 习题 5：点源 lambda=1um、2 km、r0=2 cm、望远镜 2 m
% 习题 6：点源 lambda=1um、75 km、r0=10 cm、望远镜 1 m
for cas = 1:2
    if cas == 1, Dz5 = 2e3;   r05 = 0.02; D25 = 2; nm = '习题 5';
    else,        Dz5 = 75e3;  r05 = 0.10; D25 = 1; nm = '习题 6'; end
    lam = 1e-6;
    DROI = 4*D25;  D15 = lam*Dz5/DROI;  delta1v = 2*D15/5;   % 中央瓣 5 点采样
    % 真空
    N_vac = (D15/(2*delta1v) + D25/(2*delta1v) + lam*Dz5/(2*delta1v*delta1v));
    dzmax_v = min(delta1v,delta1v)^2*2^ceil(log2(N_vac))/lam;
    nmin_v = ceil(Dz5/dzmax_v) + 1;
    % 湍流（c=2）
    D1p = D15 + 2*lam*Dz5/r05;  D2p = D25 + 2*lam*Dz5/r05;
    dn_max = (lam*Dz5 - D2p*delta1v)/D1p;                    % 约束 1
    dnv = min(dn_max, D25/25);                               % 观测面 ~25-50 点跨口径
    N_turb = (D1p/(2*delta1v) + D2p/(2*dnv) + lam*Dz5/(2*delta1v*dnv));
    Nt = 2^ceil(log2(N_turb));
    dzmax_t = min(delta1v,dnv)^2*Nt/lam;
    nmin_t = ceil(Dz5/dzmax_t) + 1;
    fprintf('%s：真空 N=2^%d、nmin=%d 次部分传播；湍流 delta_n<=%.1f mm、N=2^%d、nmin=%d\n', ...
        nm, ceil(log2(N_vac)), nmin_v, dnv*1e3, ceil(log2(N_turb)), nmin_t);
end

%% 习题 7(a)：解析 r0 与 log 振幅方差（Cn2 = 1e-17、Dz = 100 km）
lam7 = 1e-6; k7 = 2*pi/lam7; Dz7 = 100e3; Cn27 = 1e-17;
r0pw7  = (0.423*k7^2*Cn27*Dz7)^(-3/5);
r0sw7  = (0.423*k7^2*Cn27*(3/8)*Dz7)^(-3/5);
B7 = exp(gammaln(11/6)*2 - gammaln(11/3));
sch_pw = 0.563*k7^(7/6)*Cn27*Dz7^(11/6)*(6/11);      % 式 (9.63) 常数 Cn2
sch_sw = 0.563*k7^(7/6)*Cn27*Dz7^(11/6)*B7;          % 式 (9.64) 常数 Cn2
fprintf('习题 7(a)：r0,pw = %.4f m、r0,sw = %.4f m；sigma_pw^2 = %.4f、sigma_sw^2 = %.4f\n', ...
    r0pw7, r0sw7, sch_pw, sch_sw);

%% 习题 7(b)：3 块屏 3 方程（alpha = 0, 1/2, 1）的唯一解
a3 = [0 0.5 1];
A7 = [a3.^(5/3); (1-a3).^(5/6).*a3.^(5/6); (1-a3).^(5/6)];
b7 = [r0sw7^(-5/3); sch_sw/1.33*(k7/Dz7)^(5/6); sch_pw/1.33*(k7/Dz7)^(5/6)];
x7 = A7 \ b7;
fprintf('习题 7(b)：X = r0i^(-5/3) = [%.3f %.3f %.3f] -> r0i = [%.3f %.3f %.3f] m\n', ...
    x7, x7.^(-3/5));
fprintf('          回代 A*X - b 的最大偏差 = %.2e；sigma_pw^2 匹配值 = %.3f（>0.25，弱起伏界）\n', ...
    max(abs(A7*x7 - b7)), sch_pw);

%% 习题 7(c)：7 块屏的约束优化（有界最小二乘，Python 用 lsq_linear）
a7 = (0:6)/6;
A7c = [a7.^(5/3); (1-a7).^(5/6).*a7.^(5/6)];
b7c = [r0sw7^(-5/3); sch_sw/1.33*(k7/Dz7)^(5/6)];
x0 = ones(7,1)*mean(b7c);
lb = zeros(7,1);
rmax = 0.1;
ub = rmax/1.33*(k7/Dz7)^(5/6) ./ max(A7c(2,:), 1e-12);
ub(A7c(2,:) < 1e-12) = 50^(-5/3);
fun = @(X) sum((A7c*X(:) - b7c).^2);
opts = optimoptions('fmincon', 'Display', 'off');
[X7, ~] = fmincon(fun, x0, [], [], [], [], lb, ub, [], opts);
fprintf('习题 7(c)：7 屏 r0i = [%.3f %.3f %.3f %.3f %.3f %.3f %.3f] m；', X7.^(-3/5));
bp = A7c*X7;
fprintf('回代 r0,sw = %.4f m（目标 %.4f）、sigma_sw^2 = %.4f（目标 %.4f）\n', ...
    bp(1)^(-3/5), r0sw7, bp(2)*1.33*(Dz7/k7)^(5/6), sch_sw);
