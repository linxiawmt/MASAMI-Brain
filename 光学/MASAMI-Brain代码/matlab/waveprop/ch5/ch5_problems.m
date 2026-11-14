%% ===== 书 Schmidt PM199 第 5 章 =====
%  习题 1-6 复现（原书 5.3 节）：色差焦斑、像散/球差 + 离焦演化、
%  Strehl 四算法对比、环形孔径与稀疏孔径 PSF。
%  与 python/waveprop/第5章_代码实战.ipynb 工程 3 数值行为一致。
%
%  原书习题：
%    1. BK7 的 Sellmeier 色散 -> n_F/n_d/n_C -> 平凸透镜焦距与衍射受限
%       焦斑直径；再用相干成像例（5.2.1）算 PSF 切片，加多档离焦。
%       N=1024, delta=0.199 mm。
%    2. 1 波 Zernike 初级像散 + 多档离焦的 PSF 演化。
%       L=4 m, D=2 m, N=512, lambda=1 um, f=16 m。
%    3. 1 波 Zernike 初级球差 + 多档离焦的 PSF 演化（参数同习题 2）。
%    4. W = 0.07 Z4 + 0.05 Z5 - 0.05 Z6 + 0.03 Z7 - 0.03 Z8（式 5.60）的
%       Strehl：(a) 式 (5.26)+(5.54)/(5.55)；(b) 仿真 PSF（类比 5.2.1）。
%       L=8 m, D=2 m, N=512, lambda=1 um, f=64 m。
%    5. 内外径 1 m / 2 m 环形孔径 vs 满孔径 2 m 的 PSF 与 Strehl（峰值比），
%       参数同习题 4，解析确认。
%    6. 三个直径 1 m 圆孔（中心按原书印刷 (0.6,0.6), (-0.6,0.6), (0,0.6) m，
%       三孔圆心共线于 y=0.6 且两两间距 0.6 m < 直径 1 m——疑为排印笔误，
%       按原文实现）组成的稀疏孔径 PSF，解析确认（单孔谱 x 干涉因子）。
%
%  说明：习题 2-6 的 PSF 用 |ft2(P, delta)|^2（平面波 -> 焦点，即点物的
%  相干像强度），与 5.2.1 的 ASF 同一计算。MATLAB/Python 差异：数组 1 起
%  （剖面在 (N/2+1,:)）、besselj(1,.) 对应 scipy j1、zernike/zrf 与
%  ch5_zernike.m 相同（本文件复制一份以便独立运行）。
clear; clc; close all;

%% 习题 1：BK7 色散、焦距与衍射受限焦斑
lam = struct('F', 0.48612, 'd', 0.58756, 'C', 0.65627);   % um
ref = struct('F', 1.52238, 'd', 1.51680, 'C', 1.51432);   % Schott 手册六位值
R = 51.68e-3;  Dlens = 12.7e-3;           % 平凸透镜（习题 1b）
fprintf('===== 习题 1：BK7 色散、焦距与衍射受限焦斑 =====\n');
keys = {'F', 'd', 'C'};
for k = 1:3
    w = lam.(keys{k});
    n = sellmeier_bk7(w);
    f = R / (n - 1);                      % 薄平凸透镜 f = R/(n-1)
    airy = 2.44 * (w*1e-6) * f / Dlens;   % Airy 直径（第一暗环）
    f_of.(keys{k}) = f;  spot_of.(keys{k}) = airy;
    fprintf('    n(%s) = %.6f（手册 %.5f，差 %.0e）  f = %7.3f mm  Airy 直径 = %5.2f um\n', ...
        keys{k}, n, ref.(keys{k}), abs(n - ref.(keys{k})), f*1e3, airy*1e6);
end
Vd = (ref.d - 1) / (ref.F - ref.C);
fprintf('    阿贝数 V_d = (n_d-1)/(n_F-n_C) = %.2f（BK7 标准 64.17）\n', Vd);

%% 习题 1c：三波长 PSF 切片 + 离焦演化（N=1024, delta=0.199mm）
N1 = 1024;  d1m = 0.199e-3;  L1 = N1 * d1m;
fprintf('===== 习题 1c：三波长 PSF 切片 + 离焦演化 =====\n');
psf_slices = struct();
peak0_d = [];
for k = 1:3
    w = lam.(keys{k});
    [x1, y1, r1, th1, delta1] = pupil_grid(L1, N1);
    P = hypot(x1, y1) <= Dlens/2;         % 无像差光瞳
    [psf, u1] = psf_of(P, delta1, w*1e-6, f_of.(keys{k}), N1);
    row = psf(N1/2+1, :);
    half = max(row) / 2;
    above = find(row >= half);
    fwhm = (above(end) - above(1)) * (u1(2) - u1(1));
    psf_slices.(keys{k}) = {row / max(row), u1, spot_of.(keys{k})};
    if k == 2, peak0_d = max(psf(:)); end
    fprintf('    lambda=%s: PSF 半高全宽 = %5.2f um vs Airy 直径 %5.2f um\n', ...
        keys{k}, fwhm*1e6, spot_of.(keys{k})*1e6);
end
% 离焦演化（lambda_d，Zernike 离焦系数单位 RMS 波）
dfocus = [0.0, 0.125, 0.25, 0.5];
defoc_slices = cell(1, numel(dfocus));
for q = 1:numel(dfocus)
    [x1, y1, r1, th1, delta1] = pupil_grid(L1, N1);
    Wq = dfocus(q) * zernike(4, 2*r1/Dlens, th1);
    [psf, u1] = psf_of(double(hypot(x1, y1) <= Dlens/2) ...
        .* exp(1i * 2*pi * Wq), delta1, lam.d*1e-6, f_of.d, N1);
    defoc_slices{q} = {psf(N1/2+1, :) / max(psf(:)), u1};
    fprintf('    离焦 %.3f 波：PSF 峰值 = 零离焦的 %.4f\n', ...
        dfocus(q), max(psf(:)) / peak0_d);
end

%% 习题 2/3：像散、球差 + 离焦的 PSF 演化（L=4m, D=2m, f=16m）
L2 = 4.0;  N2 = 512;  D2 = 2.0;  wvl2 = 1e-6;  f2 = 16.0;
dfocus2 = [-2.0, -1.0, 0.0, 1.0, 2.0];
[x2, y2, r2, th2, delta2] = pupil_grid(L2, N2);
psf_dl2 = psf_of(double(hypot(x2, y2) <= D2/2), delta2, wvl2, f2, N2);
fprintf('===== 习题 2/3：像散、球差 + 离焦的 PSF 演化 =====\n');
names = {'像散 i=6', '球差 i=11'};  modes = [6, 11];
for g = 1:2
    for q = 1:numel(dfocus2)
        Wq = 1.0 * zernike(modes(g), 2*r2/D2, th2) ...
           + dfocus2(q) * zernike(4, 2*r2/D2, th2);
        psf = psf_of(double(hypot(x2, y2) <= D2/2) .* exp(1i * 2*pi * Wq), ...
                     delta2, wvl2, f2, N2);
        psf_grids{g}{q} = psf / max(psf(:));             %#ok<SAGROW>
        fprintf('    %s + 离焦 %+.1f 波：Strehl（PSF 峰值比）= %.4f\n', ...
            names{g}, dfocus2(q), max(psf(:)) / max(psf_dl2(:)));
    end
end

%% 习题 4：五阶混合像差的 Strehl（L=8m, D=2m, f=64m）
L4 = 8.0;  N4 = 512;  D4 = 2.0;  wvl4 = 1e-6;  f4 = 64.0;
[x4, y4, r4, th4, delta4] = pupil_grid(L4, N4);
coef = [4, 0.07; 5, 0.05; 6, -0.05; 7, 0.03; 8, -0.03];   % 式 (5.60)
W4 = zeros(N4);
for q = 1:size(coef, 1)
    W4 = W4 + coef(q, 2) * zernike(coef(q, 1), 2*r4/D4, th4);
end
sig2 = sum(coef(:, 2).^2);                % 式 (5.26)
sig2_phi = (2*pi)^2 * sig2;
ap4 = hypot(x4, y4) <= D4/2;
P_ab = ap4 .* exp(1i * 2*pi * W4);
P_dl = double(ap4);
S_pupil = abs(sum(P_ab, 'all'))^2 / abs(sum(P_dl, 'all'))^2;   % 式 (5.46)
[psf_dl4, u4] = psf_of(P_dl, delta4, wvl4, f4, N4);
psf_ab4 = psf_of(P_ab, delta4, wvl4, f4, N4);
S_sim = max(psf_ab4(:)) / max(psf_dl4(:));
fprintf('===== 习题 4：五阶混合像差的 Strehl =====\n');
fprintf('    sigma^2 = sum a_i^2 = %.4f 波^2 -> sigma_phi^2 = %.4f rad^2\n', sig2, sig2_phi);
fprintf('    (a) 式 (5.26)+(5.54)：S = 1 - sigma_phi^2 = %.4f\n', 1 - sig2_phi);
fprintf('        式 (5.55) 经验：S = exp(-sigma_phi^2) = %.4f\n', exp(-sig2_phi));
fprintf('    (b) 光瞳积分（式 5.46）：S = %.4f；仿真 PSF 峰值比：S = %.4f\n', ...
    S_pupil, S_sim);

%% 习题 5：环形孔径（内外径 1m/2m，eps=0.5）vs 满孔径
ap_ann = double(hypot(x4, y4) <= D4/2) - double(hypot(x4, y4) <= D4/4);
psf_ann = psf_of(ap_ann, delta4, wvl4, f4, N4);
S_ann = max(psf_ann(:)) / max(psf_dl4(:));
S_analytic = (1 - 0.5^2)^2;               % 轴上场 ∝ 面积：(1-eps^2)^2
c = N4/2 + 1;
row_ann = psf_ann(c, :);  row_dl = psf_dl4(c, :);
idx0  = c + find(row_ann(c:end) < row_ann(c)*1e-6, 1) - 1;
idx0_dl = c + find(row_dl(c:end) < row_dl(c)*1e-6, 1) - 1;
fprintf('===== 习题 5：环形孔径 =====\n');
fprintf('    峰值比 Strehl = %.4f vs 解析 (1-eps^2)^2 = %.4f（相对差 %.1e）\n', ...
    S_ann, S_analytic, abs(S_ann - S_analytic)/S_analytic);
fprintf('    第一暗环半径：环孔 %.2f um vs 满孔 %.2f um（环孔主瓣更窄）\n', ...
    u4(idx0)*1e6, u4(idx0_dl)*1e6);

%% 习题 6：稀疏孔径（三个 1m 圆孔，中心按原书印刷）
centers = [0.6, 0.6; -0.6, 0.6; 0.0, 0.6];   % 原书印刷：三孔共线 y=0.6
fprintf('===== 习题 6：稀疏孔径 =====\n');
% 三孔两两间距 0.6 m < 直径 1 m：按"求和"孔径做解析对照（FT 线性），
% 按"并集"孔径（物理透过率）算峰值比；两者在重叠区相差二阶小量
ap_sum = zeros(N4);  ap_sp = zeros(N4);
for q = 1:3
    dsk = double(hypot(x4 - centers(q,1), y4 - centers(q,2)) <= 0.5);
    ap_sum = ap_sum + dsk;
    ap_sp = max(ap_sp, dsk);
end
psf_sp  = psf_of(ap_sp, delta4, wvl4, f4, N4);
psf_sum = psf_of(ap_sum, delta4, wvl4, f4, N4);
% 解析：U(f) = FT{半径 0.5 圆盘}(f) * sum_j exp(-i 2pi f.c_j)
fxn = (-N4/2 : N4/2-1) / (N4 * delta4);
[FXn, FYn] = meshgrid(fxn, fxn);
rho = hypot(FXn, FYn);
U1 = 0.5 * besselj(1, 2*pi*0.5*rho) ./ max(rho, eps);
U1(rho == 0) = pi * 0.25;                 % jinc(0, R=0.5) = 面积
inter = zeros(N4);
for q = 1:3
    inter = inter + exp(-1i * 2*pi * (FXn*centers(q,1) + FYn*centers(q,2)));
end
psf_an = abs(U1 .* inter).^2;
% 数值 PSF 带 delta^4 定标；对照用归一化形状（各自除以轴上峰值）
samples = [c c; c c+1; c+2 c+3; c+5 c+2; c+9 c];
err = 0;
for q = 1:size(samples, 1)
    ii = samples(q, 1);  jj = samples(q, 2);
    err = max(err, abs(psf_sum(ii, jj)/psf_sum(c, c) - psf_an(ii, jj)/psf_an(c, c)));
end
area_ratio = sum(ap_sp(:)) * delta4^2 / (pi * 1.0^2);
S_sp = max(psf_sp(:)) / max(psf_dl4(:));
fprintf('    PSF（求和孔径）vs 解析（单孔谱 x 干涉因子）归一化形状最大差 = %.2e\n', err);
fprintf('    三孔/满孔面积比 = %.4f；峰值比 = %.4f（= 面积比平方 %.4f）\n', ...
    area_ratio, S_sp, area_ratio^2);

%% 绘图（对应 notebook 图 5-10/11/12）
figure('Position', [40 40 1400 680]);
sgtitle('习题 1：BK7 三波长焦斑（上排，红线 = ±Airy 半径）与 d 光离焦演化（下排）');
for k = 1:3
    subplot(2, 4, k);
    dat = psf_slices.(keys{k});
    row = dat{1};  u1 = dat{2};  airy = dat{3};
    m = abs(u1) < 5*airy;
    plot(u1(m)*1e6, row(m), 'k-', 'LineWidth', 1.2);  hold on;
    xline(airy*1e6, 'r:');  xline(-airy*1e6, 'r:');
    title(sprintf('lambda=%s um：Airy 直径 %.1f um', keys{k}, airy*1e6), 'FontSize', 10);
    xlabel('u [um]');  grid on;
    if k == 1, ylabel('归一化辐照度'); end
end
subplot(2, 4, 4);  axis off;              % 上排第 4 格留白
for q = 1:numel(dfocus)
    subplot(2, 4, 4 + q);
    dat = defoc_slices{q};
    m = abs(dat{2}) < 8*spot_of.d;
    plot(dat{2}(m)*1e6, dat{1}(m), 'b-', 'LineWidth', 1.2);
    title(sprintf('离焦 %.1f 波（λ=d）', dfocus(q)), 'FontSize', 10);
    xlabel('u [um]');  grid on;
end

figure('Position', [40 60 1600 860]);
for g = 1:2
    for q = 1:numel(dfocus2)
        subplot(2, 5, (g-1)*5 + q);
        imagesc(log10(psf_grids{g}{q} + 1e-6));
        caxis([-5 0]);  colormap inferno;
        xlim([-40 40]);  ylim([-40 40]);
        title(sprintf('%s, 离焦 %+.1f 波', names{g}, dfocus2(q)), 'FontSize', 10);
        xlabel('u [um]');
        if q == 1, ylabel('v [um]'); end
    end
end
sgtitle('习题 2/3：1 波像散/球差加离焦的焦斑演化（PSF，log10 显示）');

figure('Position', [40 80 1500 480]);
subplot(1, 3, 1);
imagesc(log10(psf_ann / max(psf_ann(:)) + 1e-9));  caxis([-5 0]);  colormap inferno;
xlim([-80 80]);  ylim([-80 80]);  colorbar;
title('习题 5：环形孔径 PSF（log10）');
subplot(1, 3, 2);
imagesc(log10(psf_sp / max(psf_sp(:)) + 1e-9));  caxis([-5 0]);  colormap inferno;
xlim([-80 80]);  ylim([-80 80]);  colorbar;
title('习题 6：稀疏孔径 PSF（log10）');
subplot(1, 3, 3);
mm = abs(u4) <= 120e-6;
plot(u4(mm)*1e6, row_dl(mm)/max(row_dl), 'k-', 'LineWidth', 1.3); hold on;
plot(u4(mm)*1e6, row_ann(mm)/max(row_ann), 'b-', 'LineWidth', 1.1);
plot(u4(mm)*1e6, psf_sp(c, mm)/max(psf_sp(:)), 'r-.', 'LineWidth', 1.1);
xlabel('u [um]');  ylabel('归一化 PSF');  grid on;
title('v=0 切片对比');
legend({'满孔 2 m', '环孔 ε=0.5（自归一）', '三孔稀疏（自归一）'}, 'FontSize', 8);

%% ===== 本地函数 =====
function n = sellmeier_bk7(wvl_um)
% BK7 的 Sellmeier 色散（式 5.56，Schott 手册系数，lambda 以 um 计）
B = [1.03961212, 0.231792344, 1.01046945];
C = [0.00600069867, 0.0200179144, 103.560653];
n2 = 1.0;
for q = 1:3
    n2 = n2 + B(q) * wvl_um^2 / (wvl_um^2 - C(q));
end
n = sqrt(n2);
end

function [x, y, r, theta, delta] = pupil_grid(L, N)
% 居中网格与极坐标
delta = L / N;
[x, y] = meshgrid((-N/2 : N/2-1) * delta);
r = hypot(x, y);
theta = atan2(y, x);
end

function [psf, u1d] = psf_of(P, delta, wvl, f, N)
% PSF = |ft2(P)|^2 与像面坐标（间隔 lambda f/(N delta)）
h = ft2(P, delta);
psf = abs(h).^2;
u1d = (-N/2 : N/2-1) * wvl * f / (N * delta);
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）
G = fftshift(fft2(fftshift(g))) * delta^2;
end

function Z = zernike(i, r, theta)
% Noll 单指标 Zernike（m=0 分支已按勘误补 sqrt(n+1)）；与 ch5_zernike.m 相同
[n, m] = zidx(i);
R = zrf(n, m, r);
if m == 0
    Z = sqrt(n+1) * R;
elseif mod(i, 2) == 0
    Z = sqrt(2*(n+1)) * R .* cos(m*theta);
else
    Z = sqrt(2*(n+1)) * R .* sin(m*theta);
end
end

function R = zrf(n, m, r)
% Zernike 径向多项式（式 5.5a，阶乘用 gamma 实现）
R = zeros(size(r));
for s = 0:(n-m)/2
    num = (-1)^s * gamma(n - s + 1);
    den = gamma(s+1) * gamma((n+m)/2 - s + 1) * gamma((n-m)/2 - s + 1);
    R = R + num / den * r.^(n - 2*s);
end
end

function [n, m] = zidx(i)
% 表 5.2 的 (i,n,m) 映射（Noll 单指标）；与 ch5_zernike.m 相同
persistent rows
if isempty(rows)
    rows = zeros(0, 2);
    for n0 = 0:40
        for m0 = mod(n0, 2):2:n0
            rows(end+1, :) = [n0, m0];    %#ok<AGROW>
            if m0 > 0
                rows(end+1, :) = [n0, m0];  %#ok<AGROW>
            end
        end
    end
end
n = rows(i, 1);
m = rows(i, 2);
end
