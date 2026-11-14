%% ===== 书 Schmidt PM199 第 5 章 =====
%  相干与非相干成像：ASF/PSF 与 ATF/OTF（原书 5.2 节，Listings 5.4-5.5）
%  与 python/waveprop/第5章_代码实战.ipynb 工程 2 数值行为一致。
%
%  原书要点：
%    成像系统对光是线性系统，但"线性在哪个量上"取决于照明的相干性：
%      相干照明：场线性。U_i = h ⊗ U_o（式 5.27），振幅扩散函数（ASF）
%                h = (1/(lambda zi)) F{P}（式 5.29/5.30）；
%                传递函数 H = F{h} = lambda zi * P(-lambda zi f)（式 5.34，
%                lambda zi 因子正确：二维换元雅可比 lambda zi^2 抵消
%                1/(lambda zi) 后剩一个；H(0,0) = lambda zi，量纲 m^2）。
%                圆孔截频 f0 = D/(2 lambda zi)，像是严格带限的。
%      非相干照明：强度线性。I_i = |h|^2 ⊗ I_o（式 5.35），|h|^2 即 PSF；
%                OTF：H = F{|h|^2} / ∫∫|h|^2 dudv（式 5.38）
%                   = 光瞳自相关 / ∫∫|P|^2 dxdy（式 5.40/5.41），
%                圆孔截频 2 f0（式 5.42），无像差时 H(f) >= 0。
%    波前像差只属于系统、与照明无关，但它对像的影响跟相干性有关。
%    Strehl 比（5.2.3）：S = |∫∫ P e^{i phi} dxdy|^2 / |∫∫ P dxdy|^2（式 5.46），
%    小像差时 S ≈ 1 - sigma_phi^2（Maréchal，式 5.54）≈ exp(-sigma_phi^2)（式 5.55）。
%
%  验证（编号与 notebook 工程 2 一致）：
%   (0) myconv2 正确性：FFT 线性卷积取中 vs conv2(g,h,'same')；
%   (1) 书例复现：三缝物体 0.05 波离焦的相干像与非相干像（原书图 5.5/5.6）；
%   (2) ATF 核对：F{h} 与 lambda zi*P(-lambda zi f) 逐点一致（式 5.34），
%       相干像谱 |f|>f0 泄漏（理想像严格带限，泄漏来自中心裁切窗）；
%   (3) OTF 核对：无像差圆孔数值 OTF vs 解析式 (5.42)，两条数值路线
%       （|h|^2 的 FT vs 光瞳自相关）在各自可信频段一致；截频 2f0；H>=0；
%   (4) Strehl 四算法（0.05 波离焦）；
%   (5) Parseval：sum|h|^2 du^2 = sum|P|^2 delta^2。
%
%  MATLAB/Python 差异（对应 notebook 内注释）：
%    - myconv2 用 FFT 垫零实现；对偶数尺寸方阵它与 conv2(g,h,'same') 取中
%      一致（full 卷积长 2N-1，'same' 的居中起点 = floor((K-N)/2)+1，与
%      numpy 版 0-based (K-N)//2 取同一批元素）。
%    - np.fft.fft2(g, (K,K)) 垫零 = MATLAB fft2(g, K, K)。
%    - P(-lambda zi f) 样本 = circshift(rot90(P,2), [1 1])（flip+roll 对齐
%      偶数 N 的居中反转）。
%    - cart2pol -> hypot/atan2。
clear; clc; close all;

%% ---- 原书 Listing 5.4/5.5 的参数 ----
N = 256;  L = 0.1;  D = 0.07;
delta = L / N;  wvl = 1e-6;  z = 0.25;
lam_z = wvl * z;                          % lambda zi（原书记号）
[x, y] = meshgrid((-N/2 : N/2-1) * delta);
r = hypot(x, y);  theta = atan2(y, x);
apM = hypot(x, y) <= D/2;                 % circ(x,y,D) 半径参数版

%% 验证 0：myconv2 正确性
rng(1);
gs_ = randn(6);  hs_ = randn(3);
ref = conv2(gs_, hs_, 'same');
got = real(myconv2(gs_, hs_));
fprintf('===== 验证 0：myconv2 vs conv2(same) =====\n');
fprintf('    max|myconv2 - conv2(same)| = %.2e\n', max(abs(got - ref), [], 'all'));

%% 验证 1 + 4：书例复现（0.05 波离焦）与 Strehl
W = 0.05 * zernike(4, 2*r/D, theta);      % 0.05 波 Zernike 离焦
P_ab = apM .* exp(1i * 2*pi * W);         % 广义光瞳（式 5.1）
P_dl = apM;                               % 无像差光瞳
h_ab = ft2(P_ab, delta) / lam_z;          % ASF（式 5.29，补上 1/(lambda zi)）
h_dl = ft2(P_dl, delta) / lam_z;
delta_u = lam_z / (N*delta);              % 像面采样 = lambda zi/(N delta)
u1d = (-N/2 : N/2-1) * delta_u;
[u, v] = meshgrid(u1d, u1d);
f0 = D / (2 * lam_z);                     % 相干截频
fprintf('===== 验证 1 + 4：书例复现（0.05 波离焦）与 Strehl =====\n');
fprintf('    delta_u = %.3f um；f0 = %.0f /m\n', delta_u*1e6, f0);
% 物体（Listing 5.4 第 22-23 行）：三条 50 um x 200 um 竖缝
obj = (rect((u-1.4e-4)/5e-5) + rect(u/5e-5) + rect((u+1.4e-4)/5e-5)) ...
    .* rect(v/2e-4);
img_coh = myconv2(obj, h_ab);                       % 式 (5.27)：场卷场
img_inc = myconv2(abs(obj).^2, abs(h_ab).^2);       % 式 (5.35)：强度卷强度
I_coh = abs(img_coh).^2;
I_inc = real(img_inc);
% Strehl（验证 4）：四种算法
S_pupil = abs(sum(P_ab, 'all'))^2 / abs(sum(P_dl, 'all'))^2;   % 式 (5.46) 离散
S_psf   = abs(h_ab(N/2+1, N/2+1))^2 / abs(h_dl(N/2+1, N/2+1))^2;
sig2_w  = mean(W(apM).^2);                % sigma^2（波^2）：Z4 单位 RMS -> 0.05^2
sig2_phi = (2*pi)^2 * sig2_w;
fprintf('    0.05 波 Z4 离焦：sigma_W^2 = %.4f 波^2 -> sigma_phi^2 = %.4f rad^2\n', ...
    sig2_w, sig2_phi);
fprintf('    Strehl：PSF 峰值比 = %.4f = 光瞳积分（式 5.46）%.4f\n', S_psf, S_pupil);
fprintf('            Marechal 1-sigma_phi^2（式 5.54）= %.4f\n', 1 - sig2_phi);
fprintf('            经验 exp(-sigma_phi^2)（式 5.55）= %.4f\n', exp(-sig2_phi));

%% 验证 2：ATF（式 5.34）+ 相干像带限
fx1d = (-N/2 : N/2-1) / (N * delta_u);    % 像面原生频率轴（lambda zi*f = x）
[FX, ~] = meshgrid(fx1d, fx1d);
H_num = ft2(h_ab, delta_u);
P_roll = circshift(rot90(P_ab, 2), [1 1]);   % P(-lambda zi f) 样本
dATF = max(abs(H_num - lam_z * P_roll), [], 'all');
fprintf('===== 验证 2：ATF（式 5.34）+ 相干像带限 =====\n');
fprintf('    max|F{h} - lambda zi*P(-lambda zi f)| = %.2e（相对峰值 %.1e）\n', ...
    dATF, dATF / max(abs(H_num), [], 'all'));
fprintf('    DC 值 H(0) = lambda zi = %.2e m^2（相干 ATF 不归一；OTF 才归一到 1）\n', ...
    abs(H_num(N/2+1, N/2+1)));
F_img = ft2(img_coh, delta_u);
beyond = hypot(FX, FY) > f0 * 1.001;
fprintf('    相干像谱 |f|>f0 泄漏/峰 = %.1e（理想像严格带限，泄漏来自中心裁切窗）\n', ...
    max(abs(F_img(beyond)), [], 'all') / max(abs(F_img), [], 'all'));

%% 验证 3：OTF（式 5.38/5.41/5.42）
otfA = real(ft2(abs(h_dl).^2, delta_u)) / (sum(abs(h_dl).^2, 'all') * delta_u^2);
Fp = fft2(P_dl, 2*N, 2*N);                % 垫零到 2N：位移轴可达 ±L
ac = fftshift(real(ifft2(abs(Fp).^2))) * delta^2;
xs = (-N : N-1) * delta;                  % 位移轴：-L .. L-delta（居中）
fx2 = xs / lam_z;
otfB = ac(N+1, :) / (sum(abs(P_dl).^2, 'all') * delta^2);   % fy=0 切片
% 解析式 (5.42)
otf_ana = @(f) (f <= 2*f0) .* (2/pi) .* ...
    (acos(min(max(f./(2*f0), 0), 1)) - (f./(2*f0)) .* sqrt(max(1 - (f./(2*f0)).^2, 0)));
% 路线 A 的可信频段：|h|^2 的 DFT 是循环自相关，|f| > (L-D)/(lambda zi) 即回绕
f_clean = (L - D) / lam_z;
fa  = linspace(0, 1.99*2*f0, 400);
faA = linspace(0, 0.99*f_clean, 200);
posF = fx1d >= 0;
interpA = interp1(fx1d(posF), otfA(N/2+1, posF), faA, 'linear', 'extrap');
interpB = interp1(fx2(fx2 >= 0), otfB(fx2 >= 0), fa, 'linear', 'extrap');
ana = otf_ana(fa);
anaA = otf_ana(faA);
fprintf('===== 验证 3：OTF（式 5.38/5.41/5.42） =====\n');
fprintf('    无像差圆孔 fy=0：路线A vs 解析(5.42) 干净频段（|f|<=%.2e）最大差 %.2e\n', ...
    f_clean, max(abs(interpA - anaA)));
fprintf('        路线A 超出该频段即循环回绕（对照最大差 %.2e）；路线B（垫零线性自相关）全程最大差 %.2e\n', ...
    max(abs(interp1(fx1d(posF), otfA(N/2+1, posF), fa, 'linear', 'extrap') - ana)), ...
    max(abs(interpB - ana)));
fprintf('    DC 归一：otfA(0)=%.6f，otfB(0)=%.6f（式 5.38 分母 = h 的能量）\n', ...
    otfA(N/2+1, N/2+1), otfB(N+1));
fprintf('    路线B 在 f=2f0 处 H=%.2e（截频确认）；全频段最小值 %.2e（H>=0）\n', ...
    interp1(fx2(fx2 >= 0), otfB(fx2 >= 0), 2*f0), min(otfB));
% 离焦系统的 OTF（0.05 波）——原书图 5.7 用 0.5 波展示衰减，思路相同
otf_ab = real(ft2(abs(h_ab).^2, delta_u)) / (sum(abs(h_ab).^2, 'all') * delta_u^2);

%% 验证 5：Parseval（h 与 P 能量守恒）
Eh = sum(abs(h_ab).^2, 'all') * delta_u^2;
Ep = sum(abs(P_ab).^2, 'all') * delta^2;
fprintf('===== 验证 5：Parseval =====\n');
fprintf('    sum|h|^2 du^2 = %.6e vs sum|P|^2 delta^2 = %.6e（相对差 %.2e）\n', ...
    Eh, Ep, abs(Eh - Ep)/Ep);

%% 绘图：图 5-9（自绘：像对比 + OTF 曲线 + 像谱带限）
sl = N/2 + (-60 : 59);                    % ±150 um 窗口（1-based）
figure('Position', [40 40 1500 440]);
subplot(1, 3, 1);
plot(u(1, sl)*1e3, obj(N/2+1, sl), 'k--', 'LineWidth', 1.2); hold on;
plot(u(1, sl)*1e3, I_coh(N/2+1, sl)/max(I_coh, [], 'all')*max(obj, [], 'all'), 'b-');
plot(u(1, sl)*1e3, I_inc(N/2+1, sl)/max(I_inc, [], 'all')*max(obj, [], 'all'), 'r-.');
xlabel('u [mm]');  ylabel('归一化辐照度');  grid on;
title('(a) 三缝物体的像（0.05 波离焦）');
legend({'物体辐照度', '相干像辐照度（归一）', '非相干像辐照度（归一）'}, 'FontSize', 8);
subplot(1, 3, 2);
plot(fa/(2*f0), ana, 'k-', 'LineWidth', 1.8); hold on;
plot(fa/(2*f0), interpB, 'b.--', 'MarkerSize', 4);
plot(faA/(2*f0), interpA, 'g.:', 'MarkerSize', 4);
posO = fx1d >= 0;
plot(fx1d(posO)/(2*f0), otf_ab(N/2+1, posO), 'r-.');
xline(f0/(2*f0), ':', 'Color', [0.6 0.6 0.6]);
xlabel('归一化空间频率 f/(2f0)');  ylabel('H(f)');  grid on;
title('(b) 圆孔 OTF：解析 vs 数值');
legend({'解析 (5.42) 无像差', '数值 OTF（光瞳自相关，式 5.41）', ...
    '数值 OTF（|h|^2 的FT，干净频段）', '数值 OTF（0.05 波离焦）'}, 'FontSize', 8);
ylim([-0.05 1.05]);
subplot(1, 3, 3);
F_obj = ft2(obj, delta_u);
semilogy(abs(fx1d)/f0, abs(F_obj(N/2+1, :))/max(abs(F_obj), [], 'all'), 'k-'); hold on;
semilogy(abs(fx1d)/f0, abs(F_img(N/2+1, :))/max(abs(F_img), [], 'all'), 'b-');
xline(1.0, 'r:');
xlabel('f / f0');  ylabel('归一化谱幅');  grid on;
title('(c) 相干像严格带限（f0 处截止）');
legend({'物体谱', '相干像谱'}, 'FontSize', 8, 'Location', 'southwest');
ylim([1e-8 2]);

%% ===== 本地函数 =====
function C = myconv2(g, h)
% myconv2(g, h, 1) 的语义（原书 3.1 节，该章 PDF 缺失）：全尺寸线性卷积后
% 取与 g 同尺寸的中心部分。FFT 实现：垫零到 (Ng+Nh-1) 避免循环卷积回绕，
% 再裁中间。偶数尺寸对齐：0-based 起点 (K-Ng)//2，1-based 即 (K-Ng)//2+1，
% 与 conv2(g,h,'same') 对偶数尺寸方阵取中一致。
Ng = size(g, 1);
K = Ng + size(h, 1) - 1;
full = ifft2(fft2(g, K, K) .* fft2(h, K, K));
s0 = (K - Ng)/2 + 1;                      % 1-based 裁切起点
C = full(s0:s0+Ng-1, s0:s0+Ng-1);
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）
G = fftshift(fft2(fftshift(g))) * delta^2;
end

function s = rect(x)
% rect(x)：|x|<1/2 取 1，|x|=1/2 取 1/2，其余 0
s = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
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
