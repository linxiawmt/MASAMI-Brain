%% ===== 书 Schmidt PM199 第 8 章（工程 2）：习题 1-4 复现 =====
%  原书 8.5 节习题；与 python/waveprop/ch8/ch8_problems.py 数值行为一致。
%
%  习题清单与做法：
%    1. g(x)=exp(i*pi*a^2*x^2)（式 8.26，a=4、N=128、L=4 m）的解析 FT 与
%       离散 DFT 对比；预乘 super-Gaussian 吸收边界（n=16、sigma=0.25L）后
%       重算 DFT，画实部/虚部与解析 FT 对比。
%       解析：FT{exp(i*pi*a^2*x^2)} = e^{i*pi/4}/a * exp(-i*pi*f^2/a^2)。
%       Python 实跑：低频段 |f|<=8 与解析的最大偏差：无窗 0.250（= 振幅尺度
%       1/a，混叠完全污染），加窗 0.001 —— 吸收边界把混叠贡献压掉两个量级。
%    2. 补全式 (8.20)->(8.22) 的代数步骤（alpha_2 消去），证明约束 3 与分段数
%       无关（推导打印；数值复核在工程 1 验证 0）。
%    3. 点源（lambda=1um）100 km -> 2 m 望远镜的采样图。取 delta1=delta2=1cm、
%       N=1024（第 7 章习题 2 结论，约束 4 主导）；Delta_z_max = 1.02e5 m >=
%       1e5 m（余量 2%）—— 单次传播可行，分段在这里没有收益；分段的收益
%       出现在 delta 为几十 um 的场景（书例 Delta_z_max=0.57 m << 2 m）。
%    4. 环形孔径（D_out=1.5m、D_in=0.5m、lambda=1.3um）平面波聚焦到 100 km
%       处目标平面：(a) 采样分析 + 分段数确定；(b) yn=0 切片辐照度与相位
%       （仿真 vs 解析同图）。
%       解析：源带聚焦相位 exp(-i*k*r1^2/(2*Dz)) 时菲涅耳积分的二次相位恰好
%       相消，观察面场 = e^{ikDz} e^{ikr^2/(2Dz)}/(i*lambda*Dz) x P~(f)，
%       f = r/(lambda*Dz)，P~ 为环孔光瞳的 FT：
%         P~(f) = (pi/4)*[D_out^2*2J1(pi*D_out*f)/(pi*D_out*f)
%                        - D_in^2 *2J1(pi*D_in*f) /(pi*D_in*f)]。
%       Python 实跑：delta1=50mm、delta_n=2.5mm（Airy 半径
%       1.22*lambda*Dz/D_out = 106 mm，约 40 点/半径）、N=1024、21 段
%       （每段 4.8 km <= Delta_z_max=4.9 km）；仿真/解析峰值比 0.9992，
%       ROI 内 |dI|/pk=1.0e-2、相位残差 3.6e-4 rad；
%       二次相位系数拟合 25.3 vs 解析 24.2 /m^2（差 4.7%，对应
%       拟合区间 ±8mm 内 <1e-4 rad 的相位差——系数量化对小相位差极敏感）。
%
%  MATLAB/Python 差异：MATLAB 用 besselj() 求 J1（Python 用 scipy.special.j1）；
%  索引 1-based；annulus_ft 的 f=0 极限分支要用 |u| 判断（负 f 同样走正规分支）。
clear; clc; close all;

%% ================= 习题 1：啁啾信号的 DFT 与吸收边界预乘 =================
fprintf('===== 习题 1：g(x)=exp(i*pi*a^2*x^2) 的 DFT 与吸收边界预乘 =====\n');
a = 4;  N = 128;  L = 4;
delta = L/N;
n = -N/2 : N/2-1;
x = n * delta;
g = exp(1i*pi*a^2*x.^2);                         % 式 (8.26)
% super-Gaussian 吸收边界（n=16、sigma=0.25L；r 取索引半径 |x|/delta）
sig_idx = 0.25*L/delta;
w = exp(-(abs(x)/delta).^16 / sig_idx^16);
g_win = w .* g;
% 解析 FT 与两个 DFT
f = n / (N*delta);
G_an = exp(1i*pi/4)/a * exp(-1i*pi*f.^2/a^2);
G_num = ft1(g, delta);
G_win = ft1(g_win, delta);
band = abs(f) <= 8;
fprintf('    delta=%.4f m, Nyquist=%.1f /m；混叠边界 |x| = %.2f m\n', ...
    delta, 1/(2*delta), 1/(2*a^2*delta));
fprintf('    窗在 |x|=1 m 处 |w| = %.3f；|x|=2 m 处 %.2e\n', ...
    interp1(x(x>0), w(x>0), 1.0), interp1(x(x>0), w(x>0), 2.0));
fprintf('    低频段 |f|<=8 与解析 FT 最大偏差：无窗 %.3f，加窗 %.3f（尺度 1/a=%.2f）\n', ...
    max(abs(G_num(band) - G_an(band))), max(abs(G_win(band) - G_an(band))), 1/a);
figure('Position', [60 60 1350 420]);
subplot(1, 3, 1);
plot(x, w, 'k-', 'LineWidth', 1.2); grid on;
xlabel('x [m]'); ylabel('g_{sg}(x)');
title('(a) super-Gaussian 窗（n=16, \sigma=0.25L）');
subplot(1, 3, 2);
plot(f, real(G_an), 'k-', f, real(G_num), 'r.', f, real(G_win), 'b.'); grid on;
xlim([-16 16]);
xlabel('f [m^{-1}]');
legend('解析 Re', '无窗 DFT Re', '加窗 DFT Re');
title('(b) 实部：无窗 DFT 混叠污染');
subplot(1, 3, 3);
plot(f, imag(G_an), 'k-', f, imag(G_num), 'r.', f, imag(G_win), 'b.'); grid on;
xlim([-16 16]);
xlabel('f [m^{-1}]');
legend('解析 Im', '无窗 DFT Im', '加窗 DFT Im');
title('(c) 虚部');

%% ================= 习题 2：式 (8.20)->(8.22) 的代数填充 =================
fprintf('===== 习题 2：式 (8.20) -> (8.22) 的代数填充 =====\n');
fprintf('    起点（式 8.20）：(1+Dz1/R)*d1 - w*Dz1/D1 <= d2 <= (1+Dz1/R)*d1 + w*Dz1/D1\n');
fprintf('    代入 d2 = (1-a2)*d1 + a2*dn、Dz1 = a2*Dz（式 8.21），移项合并后：\n');
fprintf('      a2*dn >= a2*d1 + a2*Dz*d1/R - w*a2*Dz/D1（左）；<= ... + w*a2*Dz/D1（右）\n');
fprintf('    两边除以 a2（>0）—— alpha_2 严格消去，得式 (8.22) == 式 (7.53)。\n');
fprintf('    数值复核见工程 1 验证 0（1000 组随机参数判定不一致 = 0 组）。\n');

%% ================= 习题 3：点源 100 km 的采样图 =================
fprintf('===== 习题 3：点源 100 km -> 2 m 望远镜的采样图 =====\n');
wvl = 1e-6;  Dz3 = 100e3;  D2p = 2;  d1p = 1e-2;  d2p = 1e-2;  Np = 1024;
fprintf('    取 delta1 = delta2 = 1 cm、N = %d（约束 4 主导）\n', Np);
fprintf('    约束 2（D1->0）：N >= %.0f；约束 4（单次）：N >= %.0f\n', ...
    D2p/(2*d2p) + wvl*Dz3/(2*d1p*d2p), wvl*Dz3/(d1p*d2p));
dzmax3 = d1p^2 * Np / wvl;                        % 式 (8.24)
fprintf('    Delta_z_max = %.3e m >= Dz = 1e5 m（余量 %.0f%%）-> 单次可行，分段无收益\n', ...
    dzmax3, (dzmax3/Dz3 - 1)*100);
fprintf('    对照书例：Delta_z_max = 0.57 m << 2 m -> 必须分段\n');
d1g = logspace(-3.3, -1.3, 300);
[d1m, dnm] = meshgrid(d1g, logspace(-3.3, -1.0, 300));
C2 = D2p./(2*dnm) + wvl*Dz3./(2*d1m.*dnm);        % 约束 2（几何项略去，D1->0）
figure('Position', [80 80 760 540]);
contourf(d1m, dnm, log2(max(C2, 1)), 6:2:18, 'LineStyle', 'none');
colormap(gray); colorbar; hold on;
plot(d1g, wvl*Dz3/Np./d1g, 'b-', 'LineWidth', 2);
plot(d1g, wvl*Dz3/Np./d1g/20, 'b--', 'LineWidth', 1.5);
plot(d1p, d2p, 'r*', 'MarkerSize', 14);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\delta_1 [m]'); ylabel('\delta_n [m]');
title('习题 3：点源 100 km 采样图（红* 上方约束 4 不紧，分段无收益）');
legend({'约束 2（log2 N）', '约束 4（单次 N=1024）', '约束 4（分 20 段）', ...
    '选点 \delta=1 cm、N=1024'}, 'Location', 'southwest');

%% ================= 习题 4：环形孔径聚焦束 100 km =================
fprintf('===== 习题 4：环形孔径（D_out=1.5m, D_in=0.5m）聚焦 100 km =====\n');
wvl4 = 1.3e-6;  Dz4 = 100e3;  D_out = 1.5;  D_in = 0.5;  D1 = D_out;
delta1 = D1/30;                                   % 30 点跨外径 -> 50 mm
deltan = 2.5e-3;                                  % Airy 半径 1.22*lambda*Dz/D = 106 mm
N4 = 1024;
D2q = 0.5;                                        % ROI 直径 0.5 m
c2 = constraint2_N(delta1, deltan, D1, D2q, wvl4, Dz4);
fprintf('    delta1=%.1f mm、delta_n=%.1f mm、约束 2：N >= %.0f -> N=%d\n', ...
    delta1*1e3, deltan*1e3, c2, N4);
constraint3_d2(delta1, D1, wvl4, Dz4, -Dz4);      % 约束 3（R=-Dz）：|1+Dz/R|=0 < D2/D1 免检
fprintf('    约束 3（R=-Dz，聚焦到目标）：|1+Dz/R|=0 < D2/D1 -> 免检\n');
dzmax4 = dz_max_partial(delta1, deltan, N4, wvl4);
nseg4 = ceil(Dz4/dzmax4);
fprintf('    约束 4：Delta_z_max = %.0f m -> 至少 %d 段（%d 平面）；取 %d 段（每段 %.1f km）\n', ...
    dzmax4, nseg4, nseg4+1, nseg4, Dz4/nseg4/1e3);
z4 = (1:nseg4) * Dz4 / nseg4;
% 源：环孔 x 聚焦相位 exp(-i k r1^2/(2 Dz))（会聚到目标面）
[x1g, y1g] = meshgrid((-N4/2 : N4/2-1) * delta1);
r1 = sqrt(x1g.^2 + y1g.^2);
ap = circle(x1g, y1g, D_out/2) - circle(x1g, y1g, D_in/2);
src = ap .* exp(-1i*pi/wvl4/Dz4 * r1.^2);
[xn4, ~, Uout4] = ang_spec_multi_prop_vac(src, wvl4, delta1, deltan, z4);
% 解析：聚焦相位与传播核啁啾相消 -> 观察面 = e^{ikDz}e^{ikr^2/2Dz}/(i w Dz) x P~(f)
xn1 = xn4(N4/2+1, :);
f4 = xn1 / (wvl4*Dz4);
Ptil = annulus_ft(f4, D_out, D_in);
U_an = exp(1i*pi/wvl4/Dz4*xn1.^2) / (1i*wvl4*Dz4) .* Ptil;
I_num = abs(Uout4(N4/2+1,:)).^2;  I_an = abs(U_an).^2;  pk = max(I_an);
win4 = abs(xn1) <= D2q/2;
eI4 = max(abs(I_num(win4) - I_an(win4)))/pk;
mph4 = win4 & (I_an > 0.05*pk);
eph4 = max(abs(angle(Uout4(N4/2+1,mph4)./U_an(mph4))));
fprintf('    解析峰值辐照度 = %.3e；仿真/解析峰值比 = %.4f\n', pk, max(I_num)/pk);
fprintf('    ROI(+-0.25 m)内：|dI|/峰值 = %.2e；相位残差(|U|>5%%峰值) = %.2e rad\n', eI4, eph4);
% 采样图 + 切片图
d1g4 = logspace(log10(delta1/5), log10(delta1*3), 260);
[d1m4, dnm4] = meshgrid(d1g4, logspace(log10(deltan/6), log10(deltan*4), 260));
C24 = constraint2_N(d1m4, dnm4, D1, D2q, wvl4, Dz4);
figure('Position', [60 60 1260 480]);
subplot(1, 2, 1);
contourf(d1m4*1e3, dnm4*1e3, log2(max(C24, 1)), 4:13, 'LineStyle', 'none');
colormap(gray); colorbar; hold on;
plot(d1g4*1e3, constraint1_d2(d1g4, D1, D2q, wvl4, Dz4)*1e3, 'k--', 'LineWidth', 1.6);
plot(delta1*1e3, deltan*1e3, 'r*', 'MarkerSize', 14);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\delta_1 [mm]'); ylabel('\delta_n [mm]');
title('(a) 习题 4(a)：采样分析（\Deltaz_{max}=4.9 km -> 21 段）');
legend({'约束 2（log2 N）', '约束 1', '选点 (50, 2.5) mm、N=1024'}, ...
    'Location', 'southwest');
subplot(1, 2, 2);
plot(xn1, I_an/pk, 'k-', xn1, I_num/pk, 'r.'); grid on;
xlim([-0.1 0.1]);
xlabel('x_n [m]'); ylabel('归一化辐照度');
legend('解析辐照度', '仿真（21 段）');
title('(b) 习题 4(b)：y_n=0 切片辐照度（主瓣 ±0.1 m）');

%% ================= 本地函数 =================
function G = ft1(g, delta)
% 原书 Listing 2.3 思路：定标正确的一维 DFT（乘 delta）。
G = fftshift(fft(fftshift(g))) * delta;
end

function out = annulus_ft(f, D_out, D_in)
% 环孔光瞳 P = circ(2r/D_out) - circ(2r/D_in) 的（径向）傅里叶变换。
% FT{circ(r/a)} = pi*a^2 * 2J1(2*pi*a*f)/(2*pi*a*f)，代入 a = D/2。
f = f(:).';
out = zeros(size(f));
for kk = 1:numel(f)
    uo = pi * D_out * f(kk);
    ui = pi * D_in * f(kk);
    if abs(uo) < 1e-12
        to = D_out^2;                            % u->0: 2J1(u)/u -> 1
    else
        to = D_out^2 * 2 * besselj(1, uo) / uo;
    end
    if abs(ui) < 1e-12
        ti = D_in^2;
    else
        ti = D_in^2 * 2 * besselj(1, ui) / ui;
    end
    out(kk) = pi/4 * (to - ti);
end
end

function [xn, yn, Uout, planes] = ang_spec_multi_prop_vac(Uin, wvl, delta1, deltan, z, sg, record)
% 任意次部分传播（书 Listing 8.1 原样整理，含吸收边界）。
if nargin < 6 || isempty(sg)
    sg = super_gaussian_absorber(size(Uin, 1));
end
if nargin < 7
    record = false;
end
N = size(Uin, 1);
[nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1));
k = 2*pi/wvl;
z = [0 z];
n = length(z);
Delta_z = z(2:n) - z(1:n-1);
alpha = z / z(n);
delta = (1-alpha) * delta1 + alpha * deltan;     % 式 (8.8)
m = delta(2:n) ./ delta(1:n-1);
x1 = nx * delta(1);
y1 = ny * delta(1);
r1sq = x1.^2 + y1.^2;
Q1 = exp(i*k/2*(1-m(1))/Delta_z(1)*r1sq);
Uin = Uin .* Q1;
planes = cell(1, n-1);
for idx = 1 : n-1
    deltaf = 1 / (N*delta(idx));
    fX = nx * deltaf;
    fY = ny * deltaf;
    fsq = fX.^2 + fY.^2;
    Z = Delta_z(idx);
    Q2 = exp(-i*pi^2*2*Z/m(idx)/k*fsq);
    Uin = sg .* ift2(Q2 .* ft2(Uin / m(idx), delta(idx)), deltaf);
    if record
        planes{idx} = Uin;
    end
end
xn = nx * delta(n);
yn = ny * delta(n);
rnsq = xn.^2 + yn.^2;
Q3 = exp(i*k/2*(m(n-1)-1)/(m(n-1)*Z)*rnsq);
Uout = Q3 .* Uin;
end

function g = super_gaussian_absorber(N, w_frac, order)
% super-Gaussian 吸收边界（式 8.1；书 Listing 8.1 第 10-12 行）。
if nargin < 2
    w_frac = 0.47;
end
if nargin < 3
    order = 16;
end
[nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1));
nsq = nx.^2 + ny.^2;
w = w_frac * N;
g = exp(-nsq.^(order/2) / w^order);
end

function [lo, hi] = constraint3_d2(d1, D1, wvl, Dz, R)
% 约束 3（式 8.22 == 式 7.53）。
slope = 1 + Dz/R;
off = wvl*Dz/D1;
lo = slope*d1 - off;
hi = slope*d1 + off;
end

function Nlo = constraint2_N(d1, dn, D1, D2, wvl, Dz)
% 约束 2（式 7.20）。
Nlo = D1./(2*d1) + D2./(2*dn) + wvl*Dz./(2*d1.*dn);
end

function dzmax = dz_max_partial(d1, dn, N, wvl)
% 式 (8.24)：Delta_z_max = min(d1,dn)^2 * N / lambda。
dzmax = min(d1, dn)^2 * N / wvl;
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）。
G = fftshift(fft2(fftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
% 原书 Listing 2.6：定标正确的二维逆 DFT（乘 (N*delta_f)^2）。
N = size(G, 1);
g = fftshift(ifft2(fftshift(G)))*(N*delta_f)^2;
end

function out = rect(x)
% rect(x)：|x|<1/2 取 1，|x|=1/2 取 1/2，其余 0。
out = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end

function out = circle(X, Y, R)
% circle(X,Y,R)：X^2+Y^2 <= R^2 内取 1，边界取 1/2（书附录 A 定义）。
r2 = X.^2 + Y.^2;
out = double(r2 < R^2) + 0.5*double(abs(r2 - R^2) < eps(R^2));
end
