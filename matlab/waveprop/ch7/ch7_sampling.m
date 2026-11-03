%% ===== 书 Schmidt PM199 第 7 章（工程 1 + 工程 2）=====
%  菲涅耳衍射的采样要求：几何带宽判据 + 二次相位判据
%  （原书 7.1-7.3 节，Listings 7.1-7.2）
%  与 python/waveprop/ch7/ch7_sampling.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准。
%
%  原书要点：
%    空间频率谱 = 平面波谱（方向余弦 alpha = lambda*f，式 7.5）；
%    传播几何限制观测孔径能看到的带宽（Coy 方法的基础）。
%    几何判据（7.2 节）：
%      alpha_max = (D1*delta2/delta1 + D2)/(2*Dz)         （式 7.12）
%      约束 1:  delta2 <= -(D2/D1)*delta1 + lambda*Dz/D1  （式 7.14）
%      约束 2:  N >= D1/(2*delta1) + D2/(2*delta2)
%                       + lambda*Dz/(2*delta1*delta2)     （式 7.20）
%    二次相位判据：f_loc = |grad phi|/(2*pi)（式 7.37）。
%    单步菲涅耳积分法：delta2 = lambda*Dz/(N*delta1)（式 7.21），
%      N >= D1*lambda*Dz/(delta1*(lambda*Dz - D2*delta1))（式 7.25 == 7.31），
%      最小距离 Dz >= D1*delta1/lambda（式 7.42，R=inf）—— 长距离有效。
%    角谱法（delta1/delta2 独立）另有两约束：
%      约束 3: (1+Dz/R)*delta1 - lambda*Dz/D1 <= delta2
%                              <= (1+Dz/R)*delta1 + lambda*Dz/D1（式 7.53）
%      约束 4: N >= lambda*Dz/(delta1*delta2)             （式 7.59）
%      —— 短距离有效。约束 4 强于约束 1+2 的组合（习题 6）。
%    一般准则：|1+Dz/R| < D2/D1 时约束 3 不起作用（式 7.60）。
%
%  原书勘误（详见笔记；本脚本按修正实现并数值验证）：
%    (a) 式 (7.46) 指数符号印反（应 +i*k/2*(...)*r1^2，与式 6.7/7.45/7.48 一致）；
%    (b) 式 (7.47) 指数符号印反（Q2[-Dz/m] 按 (6.12) 为负指数，与 Listing 6.5 一致；
%        约束 4 只用 |grad phi|，不受影响）；
%    (c) 7.3.1.2 节"发散 R<0 / 会聚 R>0"与式 (7.32) 相位定义及习题 4 矛盾：
%        R>0 发散、R<0 会聚（本脚本验证 6 数值实锤）。
%
%  Python 实跑结果（本脚本按相同参数复现）：
%    验证 0：Nmin(7.25)=65.79 -> N=128、delta2=97.7um、最小距离 8cm（书 7.43）；
%            式 7.31 与 7.25 相对差 2.2e-16；约束 2/4 的 log2 N = 8.51/8.55 -> N=512。
%    验证 1：Listing 7.1（N=128 单步）：感兴趣区域 |dI|/pk = 8.1e-3，
%            相位残差 2.1e-2 rad（原书图 7.4 复现）。
%    验证 2：Listing 7.2（N=512 角谱）：|dI|/pk = 1.2e-2，
%            相位残差 1.5e-2 rad（原书图 7.6 复现）。
%    验证 6：习题 4 会聚源 R=-0.2 m：|dI|/pk = 8.5e-3、相位残差 8.2e-3 rad；
%            输出曲率拟合 R_out=-0.107 m（几何预测 R+Dz=-0.1 m）；
%            R 符号勘误：R=+0.2 束宽 0.75->1.12 mm（发散）、R=-0.2 -> 0.38 mm（会聚）。
%
%  MATLAB/Python 差异（对应 .py 内注释）：
%    - MATLAB 1-based / numpy 0-based：居中排列的零坐标样本位置不同；
%    - meshgrid 行为一致；MATLAB 无内置 Fresnel 积分，用 integral() 数值积分
%      （Python 版用 scipy.special.fresnel）；
%    - 取 2 的幂：2^ceil(log2(Nmin))，不用 nextpow2（避免 2 的幂输入差 1）。
clear; clc; close all;

%% ================= 验证 0：约束代数总账（书例 7.43） =================
fprintf('===== 验证 0：约束代数总账（式 7.12-7.59） =====\n');
D1a = 2e-3;  D2a = 3e-3;  wvl = 1e-6;  Dza = 0.5;   % Listing 7.1 例子
d1a = D1a/50;                                       % 孔径内至少 50 点
Nmin_a = nmin_one_step(D1a, D2a, wvl, Dza, d1a);    % 式 (7.25)
N_a = 2^ceil(log2(Nmin_a));
d2a = wvl*Dza/(N_a*d1a);                            % 式 (7.21)
fprintf('    Listing 7.1 例子：D1=2mm, D2=3mm, Dz=0.5m, delta1=D1/50=%.0f um\n', d1a*1e6);
fprintf('    式 (7.25) Nmin = %.2f（书正文 66）；向上取 2 的幂 N = %d（书 128）\n', Nmin_a, N_a);
fprintf('    式 (7.21) delta2 = %.1f um（书 7.43：97.7 um）\n', d2a*1e6);
N_731 = (D1a/d1a)/(1 - D2a*d1a/(wvl*Dza));          % 式 (7.31)
fprintf('    式 (7.31)（约束 2 代入 delta2(N) 后解出）N >= %.4f，与式 (7.25) 相对差 %.1e\n', ...
    N_731, abs(N_731 - Nmin_a)/Nmin_a);
fprintf('    附加条件：lambda*Dz=%.1e > D2*delta1=%.1e（否则 N 无正解，书正文指出）\n', ...
    wvl*Dza, D2a*d1a);
fprintf('    式 (7.42) 最小距离（R=inf）= %.1f cm（书正文 8 cm）\n', ...
    min_Dz_one_step(D1a, d1a, wvl, inf)*100);
fprintf('    式 (7.41) 有限 R=+0.2 m（发散源）：Dz >= %.2f m；R=-0.2 m（会聚源）：Dz >= %.3f m\n', ...
    min_Dz_one_step(D1a, d1a, wvl, 0.2), min_Dz_one_step(D1a, d1a, wvl, -0.2));

% Listing 7.2 例子上的四条约束
D1b = 2e-3;  D2b = 4e-3;  Dzb = 0.1;
d1b = 9.4848e-6;  d2b = 28.1212e-6;
c1 = constraint1_d2(d1b, D1b, D2b, wvl, Dzb);
c2 = constraint2_N(d1b, d2b, D1b, D2b, wvl, Dzb);
[c3lo, c3hi] = constraint3_d2(d1b, D1b, wvl, Dzb, inf);
c4 = constraint4_N(d1b, d2b, wvl, Dzb);
fprintf('    Listing 7.2 例子：D1=2mm, D2=4mm, Dz=0.1m, delta1=9.4848um, delta2=28.1212um\n');
fprintf('      约束 1：delta2 <= %.1f um（选中 28.12：满足 %d）\n', c1*1e6, d2b <= c1);
fprintf('      约束 2：N >= %.2f = 2^%.2f（书正文 2^8.51）\n', c2, log2(c2));
fprintf('      约束 3：%.1f <= delta2 <= %.1f um（R=inf，选中值满足 %d）\n', ...
    c3lo*1e6, c3hi*1e6, (c3lo <= d2b && d2b <= c3hi));
fprintf('      约束 4：N >= %.1f = 2^%.2f（书正文 2^8.55）\n', c4, log2(c4));
fprintf('      N = 2^ceil(log2(max)) = %d（书 512）\n', 2^ceil(log2(max(c2, c4))));
fprintf('      式 (7.60)：|1+Dz/R|=1 < D2/D1=2 -> 约束 3 不起作用：%d\n', ...
    constraint3_is_free(D1b, D2b, wvl, Dzb, inf));

%% ================= 验证 1：Listing 7.1 复现（原书图 7.4） =================
fprintf('===== 验证 1：Listing 7.1 复现（N=128 单步传播 vs 解析） =====\n');
n = -N_a/2 : N_a/2-1;
[x1, y1] = meshgrid(n*d1a, n*d1a);
ap = rect(x1/D1a) .* rect(y1/D1a);
[x2, y2, Uout] = one_step_prop(ap, wvl, d1a, Dza);
sl = N_a/2 + 1;                                     % MATLAB 1-based：y2=0 行
x2s = x2(1, :);
Uan = fresnel_prop_square_ap(x2s, 0, D1a, wvl, Dza);
win = abs(x2s) <= D2a/2;                            % 感兴趣区域
I_num = abs(Uout(sl, :)).^2;  I_an = abs(Uan).^2;  pkA = max(I_an);
eI = max(abs(I_num(win) - I_an(win)))/pkA;
mph = win & (I_an > 0.05*pkA);
ph = angle(Uout(sl, mph)./Uan(mph));
eph = max(abs(ph));
fprintf('    N=%d, delta2=%.1f um（窗口 +-%.2f mm）\n', N_a, d2a*1e6, D2a*1e3);
fprintf('    感兴趣区域内：|dI|/峰值 = %.2e；|U|>5%%峰值处相位残差 = %.2e rad\n', eI, eph);

%% ================= 验证 2：Listing 7.2 复现（原书图 7.6） =================
fprintf('===== 验证 2：Listing 7.2 复现（N=512 角谱传播 vs 解析） =====\n');
Nb = 2^ceil(log2(max(c2, c4)));
n = -Nb/2 : Nb/2-1;
[x1, y1] = meshgrid(n*d1b, n*d1b);
apb = rect(x1/D1b) .* rect(y1/D1b);
[x2b, y2b, Uoutb] = ang_spec_prop(apb, wvl, d1b, d2b, Dzb);
slb = Nb/2 + 1;
x2s2 = x2b(1, :);
Uan2 = fresnel_prop_square_ap(x2s2, 0, D1b, wvl, Dzb);
win2 = abs(x2s2) <= D2b/2;
I_num2 = abs(Uoutb(slb, :)).^2;  I_an2 = abs(Uan2).^2;  pkB = max(I_an2);
eI2 = max(abs(I_num2(win2) - I_an2(win2)))/pkB;
mph2 = win2 & (I_an2 > 0.05*pkB);
ph2 = angle(Uoutb(slb, mph2)./Uan2(mph2));
eph2 = max(abs(ph2));
fprintf('    N=%d（约束 2=%.2f、约束 4=%.2f 的最大值向上取 2 的幂）\n', Nb, c2, c4);
fprintf('    感兴趣区域内：|dI|/峰值 = %.2e；|U|>5%%峰值处相位残差 = %.2e rad\n', eI2, eph2);

%% ================= 验证 3：图 7.5 / 图 7.7 复现（绘图） =================
fprintf('===== 验证 3：图 7.5（约束等高线）与图 7.7（式 7.60 三情形）复现 =====\n');
figure('Position', [50 50 1150 420]);
subplot(1, 2, 1);
[d1g, d2g] = meshgrid(linspace(0.2e-6, 20e-6, 200), linspace(0.2e-6, 40e-6, 200));
Nc4 = constraint4_N(d1g, d2g, wvl, Dzb);
contourf(d1g*1e6, d2g*1e6, log2(Nc4), 7:0.5:13, 'LineStyle', 'none');
colormap(gray); colorbar;
hold on;
plot(d1g(1, :)*1e6, constraint1_d2(d1g(1, :), D1b, D2b, wvl, Dzb)*1e6, 'k-.', 'LineWidth', 2);
[~, hi] = constraint3_d2(d1g(1, :), D1b, wvl, Dzb, inf);
plot(d1g(1, :)*1e6, hi*1e6, 'k--', 'LineWidth', 2);
plot(d1b*1e6, d2b*1e6, 'r*', 'MarkerSize', 12);
xlim([0 20]); ylim([0 40]);
xlabel('\delta_1 [\mum]'); ylabel('\delta_2 [\mum]');
title('(a) 约束 4 的 log2 N 等高线（对照原书图 7.5(a)）');
subplot(1, 2, 2);
Nc2 = constraint2_N(d1g, d2g, D1b, D2b, wvl, Dzb);
contourf(d1g*1e6, d2g*1e6, log2(Nc2), 7:0.5:13, 'LineStyle', 'none');
colormap(gray); colorbar;
hold on;
plot(d1g(1, :)*1e6, constraint1_d2(d1g(1, :), D1b, D2b, wvl, Dzb)*1e6, 'k-.', 'LineWidth', 2);
plot(d1g(1, :)*1e6, hi*1e6, 'k--', 'LineWidth', 2);
plot(d1b*1e6, d2b*1e6, 'r*', 'MarkerSize', 12);
xlim([0 20]); ylim([0 40]);
xlabel('\delta_1 [\mum]'); ylabel('\delta_2 [\mum]');
title('(b) 约束 2 的 log2 N 等高线（对照原书图 7.5(b)）');

figure('Position', [50 50 1150 840]);
D1c = 2e-3;  D2c = 4e-3;  Dzc = 0.1;
Rs = [Dzc/(3-1), Dzc/(0.5-1), Dzc/(-1.5-1)];        % 1+Dz/R = 3 / 0.5 / -1.5
labs = {'(a) D2/D1 < 1+\Deltaz/R', '(b) -D2/D1 < 1+\Deltaz/R < D2/D1', ...
        '(c) 1+\Deltaz/R < -D2/D1'};
d1v = linspace(0.2e-6, 20e-6, 300);
for k = 1:3
    subplot(2, 2, k);
    slope = 1 + Dzc/Rs(k);
    off = wvl*Dzc/D1c;
    c1u = constraint1_d2(d1v, D1c, D2c, wvl, Dzc);
    c3u = slope*d1v + off;
    c3l = max(slope*d1v - off, 0);
    feas = min(c1u, c3u);
    fill([d1v*1e6, fliplr(d1v*1e6)], [zeros(1, numel(d1v)), min(feas,40)*1e6], ...
        [0.6 0.6 0.6], 'EdgeColor', 'none'); hold on;
    plot(d1v*1e6, c1u*1e6, 'k-', 'LineWidth', 2);
    plot(d1v*1e6, c3u*1e6, 'k-.', 'LineWidth', 2);
    plot(d1v*1e6, c3l*1e6, 'k--', 'LineWidth', 2);
    xlim([0 20]); ylim([0 40]); grid on;
    xlabel('\delta_1 [\mum]'); ylabel('\delta_2 [\mum]');
    title(sprintf('%s：R=%.1f mm，|1+\\Deltaz/R|=%.1f vs D2/D1=2（式 7.60 判定约束 3 %s）', ...
        labs{k}, Rs(k)*1e3, abs(slope), string(constraint3_is_free(D1c, D2c, wvl, Dzc, Rs(k)))));
end
subplot(2, 2, 4); axis off;

%% ================= 验证 4：习题 1（啁啾混叠的解析证明） =================
fprintf('===== 验证 4：习题 1 —— g(x)=exp(i*pi*a^2*x^2) 的混叠解析证明 =====\n');
a1 = 4;  Nc = 128;  Lc = 4;
dc = Lc/Nc;
n = -Nc/2 : Nc/2-1;
xc = n*dc;
floc = a1^2*xc;                                     % f_loc = a^2*x（式 7.37 特例）
fnyq = 1/(2*dc);
x_alias = 1/(2*a1^2*dc);
fprintf('    delta = L/N = %.4f m；Nyquist = 1/(2*delta) = %.1f /m\n', dc, fnyq);
fprintf('    f_loc(x) = a^2*x = %.1f*x，网格边缘 |x|=%.1f m 处 f_loc = %.1f /m > Nyquist\n', ...
    a1^2, Lc/2, a1^2*Lc/2);
fprintf('    混叠边界：f_loc = Nyquist 于 |x| = %.3f m —— 恰是半网格宽\n', x_alias);
dphi = pi*a1^2*dc^2*(2*n(1:end-1) + 1);
fprintf('    相邻样本相位增量 |dphi|>pi 的样本对：%d/%d（与解析边界一致）\n', ...
    sum(abs(dphi) > pi), Nc-1);
f_app = abs(2*fnyq - floc);
fprintf('    边缘表观频率 = |2*fN - f_loc| = %.1f /m（低频回折）\n', f_app(end));

%% ================= 验证 5：习题 2/3（采样图数值） =================
fprintf('===== 验证 5：习题 2/3 —— 采样图 =====\n');
lam2 = 1e-6;  Dz2 = 100e3;  D2p = 2;
pick2 = [1e-2, 1e-2];
n2p = max(constraint4_N(pick2(1), pick2(2), lam2, Dz2), ...
          constraint2_N(pick2(1), pick2(2), 1e-6, D2p, lam2, Dz2));
fprintf('    习题 2（点源, 100 km, D2=2 m）：delta1=delta2=1 cm -> 约束 4: N>=%.0f，', ...
    constraint4_N(pick2(1), pick2(2), lam2, Dz2));
fprintf('约束 2: N>=%.0f（D1->0）-> N = 2^%d = %d\n', ...
    constraint2_N(pick2(1), pick2(2), 1e-6, D2p, lam2, Dz2), ...
    ceil(log2(n2p)), 2^ceil(log2(n2p)));
lam3 = 0.5e-6;  Dz3 = 2;  D13 = 1e-3;  D23 = 2;
pick3 = [0.4e-6, 1e-4];
n3p = max(constraint4_N(pick3(1), pick3(2), lam3, Dz3), ...
          constraint2_N(pick3(1), pick3(2), D13, D23, lam3, Dz3));
fprintf('    习题 3（D1=1 mm, 2 m, D2=2 m）：约束 1 => delta2 <= -2000*delta1+1 mm，');
fprintf('delta2>0 要求 delta1 < %.2f um；取 delta1=0.4 um, delta2=0.1 mm -> ', ...
    lam3*Dz3/D23*1e6);
fprintf('约束 4: N>=%.0f（主导），约束 2: N>=%.0f -> N = 2^%d = %d\n', ...
    constraint4_N(pick3(1), pick3(2), lam3, Dz3), ...
    constraint2_N(pick3(1), pick3(2), D13, D23, lam3, Dz3), ...
    ceil(log2(n3p)), 2^ceil(log2(n3p)));

%% ================= 验证 6：习题 4（会聚源）+ 勘误 c 实锤 =================
fprintf('===== 验证 6：习题 4 —— R=-0.2 m 会聚源：约束、仿真、解析、R 符号 =====\n');
R4 = -0.2;
rho4 = 1 + Dzb/R4;
[c3lo4, c3hi4] = constraint3_d2(d1b, D1b, wvl, Dzb, R4);
c1_4 = constraint1_d2(d1b, D1b, D2b, wvl, Dzb);
fprintf('    rho = 1+Dz/R = %.2f；约束 3：%.1f <= delta2 <= %.1f um；约束 1：delta2 <= %.1f um\n', ...
    rho4, c3lo4*1e6, c3hi4*1e6, c1_4*1e6);
fprintf('    式 (7.60)：|1+Dz/R| = %.2f < D2/D1 = %.0f -> 约束 3 不起作用：%d\n', ...
    abs(rho4), D2b/D1b, constraint3_is_free(D1b, D2b, wvl, Dzb, R4));
n = -Nb/2 : Nb/2-1;
[x1b, y1b] = meshgrid(n*d1b, n*d1b);
src4 = rect(x1b/D1b).*rect(y1b/D1b).*exp(1i*pi/wvl/R4*(x1b.^2 + y1b.^2));
[x2b4, ~, U4] = ang_spec_prop(src4, wvl, d1b, d2b, Dzb);
Uan4 = fresnel_prop_square_ap_R(x2b4(1, :), 0, D1b, wvl, Dzb, R4);
Uan_inf = fresnel_prop_square_ap_R(x2s2, 0, D1b, wvl, Dzb, inf);
fprintf('    推广解析解 R->inf 与 ch6 闭式逐位最大差 = %.1e（退化自检）\n', ...
    max(abs(Uan_inf - Uan2)));
win4 = abs(x2b4(1, :)) <= D2b/2;
I4 = abs(U4(slb, :)).^2;  Ia4 = abs(Uan4).^2;  pk4 = max(Ia4);
eI4 = max(abs(I4(win4) - Ia4(win4)))/pk4;
mph4 = win4 & (Ia4 > 0.05*pk4);
eph4 = max(abs(angle(U4(slb, mph4)./Uan4(mph4))));
fprintf('    会聚源仿真（N=%d）vs 推广解析：|dI|/峰值 = %.2e；相位残差 = %.2e rad\n', ...
    Nb, eI4, eph4);
ft = polyfit(x2b4(1, mph4).^2, unwrap(angle(U4(slb, mph4))), 1);
R_fit = pi/wvl/ft(1);
fprintf('    输出相位曲率拟合 R_out = %.4f m（几何预测 R+Dz = %.4f m）\n', R_fit, R4 + Dzb);
% 勘误 c：R=+0.2（书句称"会聚"）实际发散；R=-0.2（习题 4 称"会聚"）实际会聚
Nw = 256;  Lw = 8e-3;  dw = Lw/Nw;  w0 = 1.5e-3;
nw = -Nw/2 : Nw/2-1;
[xw, yw] = meshgrid(nw*dw, nw*dw);
for Rk = [0.2, -0.2]
    Uw = exp(-(xw.^2 + yw.^2)/w0^2).*exp(1i*pi/wvl/Rk*(xw.^2 + yw.^2));
    Iz0 = abs(Uw(Nw/2+1, :)).^2;
    w0_rms = sqrt(sum(Iz0.*xw(Nw/2+1, :).^2)/sum(Iz0));
    [~, ~, Uwp] = ang_spec_prop(Uw, wvl, dw, dw, 0.1);
    Iz1 = abs(Uwp(Nw/2+1, :)).^2;
    w1_rms = sqrt(sum(Iz1.*xw(Nw/2+1, :).^2)/sum(Iz1));
    mp = Iz1 > 0.05*max(Iz1);
    ftk = polyfit(xw(Nw/2+1, mp).^2, unwrap(angle(Uwp(Nw/2+1, mp))), 1);
    R_out = pi/wvl/ftk(1);
    fprintf('    [R=%+.1f] 高斯束 RMS 宽度 %.2f mm -> %.2f mm（%s）；输出曲率 R_out = %+.3f m\n', ...
        Rk, w0_rms*1e3, w1_rms*1e3, string(w1_rms > w0_rms), R_out);
end

%% ================= 验证 7：习题 6（约束 4 强于约束 1+2 组合） =================
fprintf('===== 验证 7：习题 6 —— 约束 4 强于约束 1+2 的组合 =====\n');
fprintf('    代数：约束 1 <=> D1*delta2 + D2*delta1 <= lambda*Dz\n');
fprintf('          两边除以 2*delta1*delta2：几何项 <= lambda*Dz/(2*delta1*delta2)\n');
fprintf('          故约束 2 的右端 <= lambda*Dz/(delta1*delta2) = 约束 4 的右端\n');
rng(7, 'twister');
viol = 0;
for it = 1:1000
    d1r = 10^(-6 + 2*rand);
    d2r = 10^(-6 + 2*rand);
    D1r = 10^(-4 + 2*rand);
    D2r = 10^(-4 + 2*rand);
    Dzr = 10^(-1 + 3*rand);
    geo = D1r/(2*d1r) + D2r/(2*d2r);
    if geo <= constraint4_N(d1r, d2r, wvl, Dzr)/2 + 1e-9
        viol = max(viol, constraint2_N(d1r, d2r, D1r, D2r, wvl, Dzr) ...
                      - constraint4_N(d1r, d2r, wvl, Dzr));
    end
end
fprintf('    1000 组随机参数（约束 1 成立域内）：max[约束2右端 - 约束4右端] = %.3e（<=0 恒成立）\n', viol);

%% ================= 绘图：图 7-8 对照（可选） =================
figure('Position', [60 60 1150 760]);
subplot(2, 2, 1);
plot(x2s*1e3, I_num/pkA, 'k-', x2s*1e3, I_an/pkA, 'r--', 'LineWidth', 1.2);
xlim([-1.5 1.5]); grid on;
xlabel('x_2 [mm]'); ylabel('归一化辐照度'); legend('数值', '解析');
title(sprintf('(a) 单步法（Listing 7.1）：N=128，|dI|/pk = %.1e', eI));
subplot(2, 2, 2);
plot(x2s*1e3, angle(Uout(sl, :)), 'k.', x2s*1e3, angle(Uan), 'r-', 'LineWidth', 1);
xlim([-1.5 1.5]); grid on;
xlabel('x_2 [mm]'); ylabel('相位 [rad]');
title(sprintf('(b) 单步法相位：残差 %.1e rad', eph));
subplot(2, 2, 3);
plot(x2s2*1e3, I_num2/pkB, 'k-', x2s2*1e3, I_an2/pkB, 'r--', 'LineWidth', 1.2);
xlim([-2.5 2.5]); grid on;
xlabel('x_2 [mm]'); ylabel('归一化辐照度'); legend('数值', '解析');
title(sprintf('(c) 角谱法（Listing 7.2）：N=512，|dI|/pk = %.1e', eI2));
subplot(2, 2, 4);
plot(x2s2*1e3, angle(Uoutb(slb, :)), 'k.', x2s2*1e3, angle(Uan2), 'r-', 'LineWidth', 1);
xlim([-2.5 2.5]); grid on;
xlabel('x_2 [mm]'); ylabel('相位 [rad]');
title(sprintf('(d) 角谱法相位：残差 %.1e rad', eph2));

%% ================= 本地函数 =================
function d2 = constraint1_d2(d1, D1, D2, wvl, Dz)
% 约束 1（式 7.14）：delta2 的上界
d2 = -D2/D1.*d1 + wvl*Dz/D1;
end

function Nlo = constraint2_N(d1, d2, D1, D2, wvl, Dz)
% 约束 2（式 7.20）：观测面网格点数下界
Nlo = D1./(2*d1) + D2./(2*d2) + wvl*Dz./(2*d1.*d2);
end

function [lo, hi] = constraint3_d2(d1, D1, wvl, Dz, R)
% 约束 3（式 7.53）：delta2 的 [下界, 上界]
slope = 1 + Dz/R;
off = wvl*Dz/D1;
lo = slope*d1 - off;
hi = slope*d1 + off;
end

function Nlo = constraint4_N(d1, d2, wvl, Dz)
% 约束 4（式 7.59）
Nlo = wvl*Dz./(d1.*d2);
end

function out = constraint3_is_free(D1, D2, wvl, Dz, R)
% 式 7.60：|1+Dz/R| < D2/D1 时约束 3 不起作用
out = abs(1 + Dz/R) < D2/D1;
end

function Nmin = nmin_one_step(D1, D2, wvl, Dz, d1)
% 单步法网格点数下界（式 7.25 == 7.31）；要求 lambda*Dz > D2*delta1
Nmin = D1*wvl*Dz/(d1*(wvl*Dz - D2*d1));
end

function Dzmin = min_Dz_one_step(D1, d1, wvl, R)
% 单步法最小传播距离（式 7.41 有限 R / 式 7.42 无限 R）
% 注意：R 按式 (7.32) 相位 exp(i*k*r1^2/(2*R)) 定义，R<0 会聚（勘误 c）
if isinf(R)
    Dzmin = D1*d1/wvl;
else
    Dzmin = D1*d1*R/(wvl*R - D1*d1);
end
end

function U = fresnel_prop_square_ap(x2, y2, D, wvl, Dz)
% 边长 D 方孔的 Fresnel 衍射解析解（式 6.1 可分离积分，省略 e^{ikDz}）：
% 观测面啁啾与配方相位相消，U = I(x2)I(y2)/(i*lambda*Dz)
s = sqrt(2/(wvl*Dz));
[tX1, tX2] = deal(s*(x2 - D/2), s*(x2 + D/2));
[tY1, tY2] = deal(s*(y2 - D/2), s*(y2 + D/2));
[Sx1, Cx1] = fresnel_int(tX1);
[Sx2, Cx2] = fresnel_int(tX2);
[Sy1, Cy1] = fresnel_int(tY1);
[Sy2, Cy2] = fresnel_int(tY2);
Ix = sqrt(wvl*Dz/2)*((Cx2 - Cx1) + 1i*(Sx2 - Sx1));
Iy = sqrt(wvl*Dz/2)*((Cy2 - Cy1) + 1i*(Sy2 - Sy1));
U = Ix.*Iy/(1i*wvl*Dz);
end

function U = fresnel_prop_square_ap_R(x2, y2, D, wvl, Dz, R)
% 方孔 + 二次相位源 exp(i*k*r1^2/(2*R)) 的解析解（习题 4a）：
%   U = exp(i*pi*(x2^2+y2^2)/(lambda*Dz)) * I(x2)*I(y2)/(i*lambda*Dz)
%   I(x) = exp(-i*pi*x^2/(rho*lambda*Dz))*sqrt(lambda*Dz/(2*rho))
%          * {[C(t2)-C(t1)] + i[S(t2)-S(t1)]}
%   rho = 1 + Dz/R，t = sqrt(2*rho/(lambda*Dz))*(x1 - x/rho)，x1 = +-D/2
% R->inf 时观测面啁啾与配方相位相消，逐位退化为 fresnel_prop_square_ap。
% 输出曲率 R_out = R + Dz（球面波几何）。
rho = 1 + Dz/R;
s = sqrt(2*rho/(wvl*Dz));
[tX1, tX2] = deal(s*(-D/2 - x2/rho), s*(D/2 - x2/rho));
[tY1, tY2] = deal(s*(-D/2 - y2/rho), s*(D/2 - y2/rho));
[Sx1, Cx1] = fresnel_int(tX1);
[Sx2, Cx2] = fresnel_int(tX2);
[Sy1, Cy1] = fresnel_int(tY1);
[Sy2, Cy2] = fresnel_int(tY2);
Ix = exp(-1i*pi*x2.^2/(rho*wvl*Dz))*sqrt(wvl*Dz/(2*rho)) ...
    *((Cx2 - Cx1) + 1i*(Sx2 - Sx1));
Iy = exp(-1i*pi*y2.^2/(rho*wvl*Dz))*sqrt(wvl*Dz/(2*rho)) ...
    *((Cy2 - Cy1) + 1i*(Sy2 - Sy1));
U = exp(1i*pi*(x2.^2 + y2.^2)/(wvl*Dz)).*Ix.*Iy/(1i*wvl*Dz);
end

function [S, C] = fresnel_int(t)
% Fresnel 积分 S(t)=int_0^t sin(pi u^2/2)du、C(t)=int_0^t cos(pi u^2/2)du。
% MATLAB 基础版无内置（符号工具箱除外）；用 integral() 数值积分
% （Python 版用 scipy.special.fresnel）。t 为数组，逐元素求值。
S = arrayfun(@(tt) integral(@(u) sin(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
C = arrayfun(@(tt) integral(@(u) cos(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
end

function [x2, y2, Uout] = one_step_prop(Uin, wvl, d1, Dz)
% 书 Listing 6.1：单步 Fresnel 传播（式 6.5；省略全局相位 e^{ikDz}，原书约定）
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
d2 = wvl*Dz/(N*d1);                            % 式 (6.16)
[x2, y2] = meshgrid(n*d2, n*d2);
Uout = 1/(1i*wvl*Dz) ...
    .* exp(1i*kk/(2*Dz)*(x2.^2 + y2.^2)) ...
    .* ft2(Uin.*exp(1i*kk/(2*Dz)*(x1.^2 + y1.^2)), d1);
end

function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
% 书 Listing 6.5：角谱传播（式 6.50 的 m 分支）：
%   Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin / m, d1), df1)
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
r1sq = x1.^2 + y1.^2;
df1 = 1/(N*d1);
[fX, fY] = meshgrid(n*df1, n*df1);
fsq = fX.^2 + fY.^2;
m = d2/d1;
[x2, y2] = meshgrid(n*d2, n*d2);
r2sq = x2.^2 + y2.^2;
Q1 = exp(1i*kk/2*(1 - m)/Dz*r1sq);
Q2 = exp(-1i*pi^2*2*Dz/m/kk*fsq);              % 式 (7.47) 勘误 b：负指数
Q3 = exp(1i*kk/2*(m - 1)/(m*Dz)*r2sq);
Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin/m, d1), df1);
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）
G = fftshift(fft2(fftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
% 原书 Listing 2.6：定标正确的二维逆 DFT（乘 (N*delta_f)^2）
N = size(G, 1);
g = fftshift(ifft2(fftshift(G)))*(N*delta_f)^2;
end

function out = rect(x)
% rect(x)：|x|<1/2 取 1，|x|=1/2 取 1/2，其余 0
out = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end
