%% ===== 书 Schmidt PM199 第 7 章（工程 3）：习题复现 =====
%  习题 1-6 复现（原书 7.4 节）
%  与 python/waveprop/ch7/ch7_sampling.py（验证 4-7）数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准。
%
%  原书习题：
%    1. g(x)=exp(i*pi*a^2*x^2)，a=4，N=128，L=4 m：不做 FT，解析证明混叠。
%    2. 点源采样图：lambda=1 um，传播 100 km，接收孔径 2 m。
%    3. 采样图：lambda=0.5 um，D1=1 mm，传播 2.0 m，接收孔径 2 m。
%    4. 把 Listing 7.2/B.5 改成会聚/发散源（式 7.62）：
%       (a) 推导带 R 的方孔解析解；(b) R=-0.2 m 画约束图；(c) 仿真 vs 解析。
%    5. 图示式 (7.60) <=> 几何光束包含在直径 D2 内（发散源 + 会聚源
%       焦点在观测面前/后四种情形的光线图）。
%    6. 代数证明约束 4 强于约束 1+2 的组合。
%
%  Python 实跑结果（本脚本按相同参数复现）：
%    习题 1：f_loc = a^2*x，Nyquist = 16/m，混叠边界 |x| = 1.000 m（半网格）；
%            相邻样本相位增量 |dphi|>pi 的样本对 63/127，边界一致。
%    习题 2：delta1=delta2=1 cm -> 约束 4: N>=1000（主导）、约束 2: N>=600
%            -> N = 2^10 = 1024。
%    习题 3：约束 1 => delta2 <= -2000*delta1 + 1 mm，delta2>0 要求 delta1
%            < 0.50 um；取 delta1=0.4 um、delta2=0.1 mm -> 约束 4: N>=25000
%            （主导）、约束 2: N>=23750 -> N = 2^15 = 32768。
%    习题 4：rho=0.5；约束 1 上界 31.0 um（主导）、约束 3 [-45.3, 54.7] um；
%            选 delta1=9.4848 um、delta2=28.1212 um、N=512：
%            |dI|/pk = 8.5e-3、相位残差 8.2e-3 rad；
%            输出曲率拟合 R_out = -0.107 m（几何预测 R+Dz = -0.1 m）；
%            勘误 c：R=+0.2 束宽 0.75->1.12 mm（发散）、R=-0.2 -> 0.38 mm（会聚）。
%    习题 6：1000 组随机参数（约束 1 成立域内）max[约束2 - 约束4] = 0。
%
%  习题 4(a) 推导（式 7.36 算符链，笔记"拓展"给出完整过程）：
%    U(x2,y2) = e^{i*pi*(x2^2+y2^2)/(lambda*Dz)} I(x2) I(y2)/(i*lambda*Dz)
%    I(x) = e^{-i*pi*x^2/(rho*lambda*Dz)} sqrt(lambda*Dz/(2*rho))
%           * {[C(t2)-C(t1)] + i[S(t2)-S(t1)]}，rho = 1 + Dz/R，
%    t = sqrt(2*rho/(lambda*Dz))*(x1 - x/rho)，x1 = +-D/2。
%    R->inf 退化为书式 (1.60) 的实现（逐位核对 4.4e-16）；
%    输出波前曲率 R_out = R + Dz（球面波几何）。
%
%  MATLAB/Python 差异：同 ch7_sampling.m 头注。
clear; clc; close all;

wvl = 1e-6;

%% ================= 习题 1：啁啾混叠的解析证明 =================
fprintf('===== 习题 1：g(x)=exp(i*pi*a^2*x^2), a=4, N=128, L=4 m =====\n');
a1 = 4;  Nc = 128;  Lc = 4;
dc = Lc/Nc;
n = -Nc/2 : Nc/2-1;
xc = n*dc;
floc = a1^2*xc;                       % 式 (7.37)：f_loc = (1/2pi) dphi/dx
fnyq = 1/(2*dc);
x_alias = 1/(2*a1^2*dc);
fprintf('    delta = L/N = %.4f m；Nyquist = 1/(2*delta) = %.1f /m\n', dc, fnyq);
fprintf('    f_loc(x) = a^2*x = %.1f*x；网格边缘 |x|=2 m 处 f_loc=%.1f /m > Nyquist\n', ...
    a1^2, a1^2*Lc/2);
fprintf('    混叠边界 |x| = %.3f m —— 恰是半网格宽，一半样本混叠\n', x_alias);
dphi = pi*a1^2*dc^2*(2*n(1:end-1) + 1);   % 相邻样本相位增量（无需 FT）
fprintf('    |dphi|>pi 的样本对：%d/%d（与解析边界一致）\n', sum(abs(dphi) > pi), Nc-1);
fprintf('    边缘表观频率 |2*fN - f_loc| = %.1f /m（低频回折）\n', abs(2*fnyq - floc(end)));
g = exp(1i*pi*a1^2*xc.^2);
figure('Position', [60 60 1050 380]);
subplot(1, 2, 1);
plot(xc, real(g), 'k-'); hold on;
xspan(x_alias, Lc/2, [1 0.8 0.8]);
xspan(-Lc/2, -x_alias, [1 0.8 0.8]);
xlabel('x [m]'); ylabel('Re g');
title('(a) 采样后 Re g：红色区域 f_{loc} > Nyquist（混叠）');
subplot(1, 2, 2);
plot(xc, floc, 'b--', xc, abs(2*fnyq - floc), 'g-.', 'LineWidth', 1.2); hold on;
yline(fnyq, 'r:');
ylim([0 70]); xlabel('x [m]'); ylabel('频率 [1/m]');
legend('f_{loc}=a^2x', '|2f_N-f_{loc}|', 'Nyquist');
title('(b) 局部频率与回折：|x|>1 m 折叠回低频');

%% ================= 习题 2/3：采样图数值 + 绘图 =================
fprintf('===== 习题 2/3：采样图 =====\n');
lam2 = 1e-6;  Dz2 = 100e3;  D2p = 2;              % 习题 2：点源
pick2 = [1e-2, 1e-2];
[d1p, d2p] = meshgrid(logspace(-4, -1, 200), logspace(-4, -1, 200));
Nc4p = constraint4_N(d1p, d2p, lam2, Dz2);
Nc2p = constraint2_N(d1p, d2p, 1e-6, D2p, lam2, Dz2);
lam3 = 0.5e-6;  Dz3 = 2;  D13 = 1e-3;  D23 = 2;   % 习题 3
pick3 = [0.4e-6, 1e-4];
[d1q, d2q] = meshgrid(logspace(-6.3, -3, 200), logspace(-5, -2, 200));
Nc4q = constraint4_N(d1q, d2q, lam3, Dz3);
Nc2q = constraint2_N(d1q, d2q, D13, D23, lam3, Dz3);
n2p = max(constraint4_N(pick2(1), pick2(2), lam2, Dz2), ...
          constraint2_N(pick2(1), pick2(2), 1e-6, D2p, lam2, Dz2));
n3p = max(constraint4_N(pick3(1), pick3(2), lam3, Dz3), ...
          constraint2_N(pick3(1), pick3(2), D13, D23, lam3, Dz3));
fprintf('    习题 2：delta1=delta2=1 cm -> 约束 4: N>=%.0f（主导），约束 2: N>=%.0f -> N=2^%d=%d\n', ...
    constraint4_N(pick2(1), pick2(2), lam2, Dz2), ...
    constraint2_N(pick2(1), pick2(2), 1e-6, D2p, lam2, Dz2), ...
    ceil(log2(n2p)), 2^ceil(log2(n2p)));
fprintf('    习题 3：约束 1 => delta1 < %.2f um；取 0.4 um/0.1 mm -> 约束 4: N>=%.0f（主导）', ...
    lam3*Dz3/D23*1e6, constraint4_N(pick3(1), pick3(2), lam3, Dz3));
fprintf('，约束 2: N>=%.0f -> N=2^%d=%d\n', ...
    constraint2_N(pick3(1), pick3(2), D13, D23, lam3, Dz3), ...
    ceil(log2(n3p)), 2^ceil(log2(n3p)));
figure('Position', [60 60 1150 460]);
subplot(1, 2, 1);
contourf(d1p, d2p, log2(max(Nc4p, 1e-6)), 9:1:14, 'LineStyle', 'none');
colormap(gray); colorbar; set(gca, 'XScale', 'log', 'YScale', 'log');
hold on; plot(pick2(1), pick2(2), 'r*', 'MarkerSize', 12);
xlabel('\delta_1 [m]'); ylabel('\delta_2 [m]');
title('(a) 习题 2：点源 100 km，D2=2 m（约束 4 的 log2 N）');
subplot(1, 2, 2);
contourf(d1q, d2q, log2(max(Nc2q, 1e-6)), 9:0.5:12, 'LineStyle', 'none');
colormap(gray); colorbar; set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
dl = linspace(0.35e-6, 1e-6, 100);
plot(dl, constraint1_d2(dl, D13, D23, lam3, Dz3), 'r-', 'LineWidth', 1.5);
plot(pick3(1), pick3(2), 'r*', 'MarkerSize', 12);
xlabel('\delta_1 [m]'); ylabel('\delta_2 [m]');
title('(b) 习题 3：D1=1 mm，2 m，D2=2 m（约束 2 等高线 + 约束 1 红线）');

%% ================= 习题 4：会聚源（式 7.62） =================
fprintf('===== 习题 4：R=-0.2 m 会聚源 =====\n');
D1b = 2e-3;  D2b = 4e-3;  Dzb = 0.1;  R4 = -0.2;
d1b = 9.4848e-6;  d2b = 28.1212e-6;
rho4 = 1 + Dzb/R4;
[c3lo4, c3hi4] = constraint3_d2(d1b, D1b, wvl, Dzb, R4);
c1_4 = constraint1_d2(d1b, D1b, D2b, wvl, Dzb);
Nmin4 = max(constraint2_N(d1b, d2b, D1b, D2b, wvl, Dzb), ...
            constraint4_N(d1b, d2b, wvl, Dzb));
N4 = 2^ceil(log2(Nmin4));
fprintf('    rho = %.2f；约束 1：delta2 <= %.1f um（主导）；约束 3：[%.1f, %.1f] um；N = %d\n', ...
    rho4, c1_4*1e6, c3lo4*1e6, c3hi4*1e6, N4);
n = -N4/2 : N4/2-1;
[x1b, y1b] = meshgrid(n*d1b, n*d1b);
src4 = rect(x1b/D1b).*rect(y1b/D1b).*exp(1i*pi/wvl/R4*(x1b.^2 + y1b.^2));
[x2b4, ~, U4] = ang_spec_prop(src4, wvl, d1b, d2b, Dzb);
Uan4 = fresnel_prop_square_ap_R(x2b4(1, :), 0, D1b, wvl, Dzb, R4);
slb = N4/2 + 1;
win4 = abs(x2b4(1, :)) <= D2b/2;
I4 = abs(U4(slb, :)).^2;  Ia4 = abs(Uan4).^2;  pk4 = max(Ia4);
eI4 = max(abs(I4(win4) - Ia4(win4)))/pk4;
mph4 = win4 & (Ia4 > 0.05*pk4);
eph4 = max(abs(angle(U4(slb, mph4)./Uan4(mph4))));
ftk = polyfit(x2b4(1, mph4).^2, unwrap(angle(U4(slb, mph4))), 1);
R_fit = pi/wvl/ftk(1);
fprintf('    仿真 vs 推广解析：|dI|/峰值 = %.2e；相位残差 = %.2e rad\n', eI4, eph4);
fprintf('    输出曲率拟合 R_out = %.4f m（几何预测 R+Dz = %.4f m）\n', R_fit, R4 + Dzb);
figure('Position', [60 60 1150 760]);
subplot(2, 2, 1);
d1v = linspace(0.2e-6, 20e-6, 300);
[lo4, hi4] = constraint3_d2(d1v, D1b, wvl, Dzb, R4);
c1v = constraint1_d2(d1v, D1b, D2b, wvl, Dzb);
fill([d1v*1e6, fliplr(d1v*1e6)], ...
     [zeros(1, numel(d1v)), min(min(c1v, hi4), 40)*1e6], [0.6 0.6 0.6], ...
     'EdgeColor', 'none'); hold on;
plot(d1v*1e6, c1v*1e6, 'k-', 'LineWidth', 2);
plot(d1v*1e6, hi4*1e6, 'k-.', 'LineWidth', 2);
plot(d1v*1e6, max(lo4, 0)*1e6, 'k--', 'LineWidth', 2);
plot(d1b*1e6, d2b*1e6, 'r*', 'MarkerSize', 12);
xlim([0 20]); ylim([0 40]); grid on;
xlabel('\delta_1 [\mum]'); ylabel('\delta_2 [\mum]');
title('(a) 习题 4(b)：R=-0.2 m 的约束 1/3（红* = 选中网格）');
subplot(2, 2, 2);
plot(x2b4(1, :)*1e3, I4/pk4, 'k-', x2b4(1, :)*1e3, Ia4/pk4, 'r--', 'LineWidth', 1.2);
xlim([-2 2]); grid on;
xlabel('x_2 [mm]'); ylabel('归一化辐照度'); legend('仿真', '解析（习题 4a）');
title(sprintf('(b) 会聚源角谱仿真 vs 推广解析：|dI|/pk = %.1e', eI4));
subplot(2, 2, 3);
plot(x2b4(1, :)*1e3, angle(U4(slb, :)), 'k.', x2b4(1, :)*1e3, angle(Uan4), 'r-', ...
    'LineWidth', 1);
xlim([-2 2]); grid on;
xlabel('x_2 [mm]'); ylabel('相位 [rad]');
title(sprintf('(c) 输出相位（曲率拟合 R_{out}=%.3f m ≈ R+\\Deltaz=%.1f m）', R_fit, R4 + Dzb));

% 习题 5（图示式 7.60）：三种 R 的光线图 —— 几何光束是否含在 D2 内
figure('Position', [60 60 1150 380]);
Rs = [Dzb/(3-1), Dzb/(0.5-1), Dzb/(-1.5-1)];             % 1+Dz/R = 3 / 0.5 / -1.5
for kk = 1:3
    subplot(1, 3, kk);
    Rk = Rs(kk);
    % 光线包络：源面 +-D1/2 出发，斜率 = (x1/R + (x2-x1)/Dz 的近轴光线)
    x1e = [-D1b/2, D1b/2];
    for e = 1:2
        % 从源面边缘出发的两条光线：发散/会聚方向
        slope = x1e(e)/Rk;                 % 球面波：光线斜率 = x1/R
        x_end = x1e(e) + slope*Dzb;
        plot([0, Dzb], [x1e(e), x_end]*1e3, 'b-'); hold on;
        plot([0, Dzb], [-x1e(e), -x_end]*1e3, 'b-');
    end
    yline(D2b/2*1e3, 'r--'); yline(-D2b/2*1e3, 'r--');
    xlim([0 Dzb]); ylim([-2.5 2.5]); grid on;
    xlabel('z [m]'); ylabel('x [mm]');
    title(sprintf('1+\\Deltaz/R=%+.1f（式 7.60：约束 3 %s）', ...
        1 + Dzb/Rk, string(constraint3_is_free(D1b, D2b, wvl, Dzb, Rk))));
end

%% ================= 习题 6：约束 4 强于约束 1+2 的组合 =================
fprintf('===== 习题 6：代数证明 + 数值核对 =====\n');
fprintf('    约束 1 <=> D1*delta2 + D2*delta1 <= lambda*Dz\n');
fprintf('    两边除以 2*delta1*delta2：D1/(2*delta1)+D2/(2*delta2) <= lambda*Dz/(2*delta1*delta2)\n');
fprintf('    故约束 2 右端 = 几何项 + lambda*Dz/(2*d1*d2) <= lambda*Dz/(d1*d2) = 约束 4 右端\n');
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
fprintf('    1000 组随机参数：max[约束2右端 - 约束4右端] = %.3e（<=0 恒成立）\n', viol);

%% ================= 本地函数 =================
function d2 = constraint1_d2(d1, D1, D2, wvl, Dz)
d2 = -D2/D1.*d1 + wvl*Dz/D1;                 % 约束 1（式 7.14）
end

function Nlo = constraint2_N(d1, d2, D1, D2, wvl, Dz)
Nlo = D1./(2*d1) + D2./(2*d2) + wvl*Dz./(2*d1.*d2);   % 约束 2（式 7.20）
end

function [lo, hi] = constraint3_d2(d1, D1, wvl, Dz, R)
slope = 1 + Dz/R;
off = wvl*Dz/D1;
lo = slope*d1 - off;
hi = slope*d1 + off;                          % 约束 3（式 7.53）
end

function Nlo = constraint4_N(d1, d2, wvl, Dz)
Nlo = wvl*Dz./(d1.*d2);                       % 约束 4（式 7.59）
end

function out = constraint3_is_free(D1, D2, wvl, Dz, R)
out = abs(1 + Dz/R) < D2/D1;                  % 式 7.60
end

function U = fresnel_prop_square_ap_R(x2, y2, D, wvl, Dz, R)
% 方孔 + 二次相位源 exp(i*k*r1^2/(2*R)) 的解析解（习题 4a）
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
S = arrayfun(@(tt) integral(@(u) sin(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
C = arrayfun(@(tt) integral(@(u) cos(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
end

function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
% 书 Listing 6.5：角谱传播（式 6.50 的 m 分支）
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
Q2 = exp(-1i*pi^2*2*Dz/m/kk*fsq);
Q3 = exp(1i*kk/2*(m - 1)/(m*Dz)*r2sq);
Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin/m, d1), df1);
end

function G = ft2(g, delta)
G = fftshift(fft2(fftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
N = size(G, 1);
g = fftshift(ifft2(fftshift(G)))*(N*delta_f)^2;
end

function out = rect(x)
out = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end

function xspan(x0, x1, c)
% 竖直阴影带（axspan 的简化版，作用于当前坐标轴）
yl = ylim;
fill([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], c, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
