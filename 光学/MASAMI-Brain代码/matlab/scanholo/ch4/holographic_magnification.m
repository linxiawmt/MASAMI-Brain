% ===== 书 Example 4.1 · 全息放大率（式 4.2-1 ~ 4.2-5） =====
% 三点物体 (4.2-1)：δ(x,y,z-z0) + δ(x-x0,y,z-z0) + δ(x,y,z-(z0+Δz0))
% 正弦全息图 (4.2-2)：
%   H_3p(x,y) ~ sin[(k0/2z0)(x²+y²)] + sin{(k0/2z0)[(x-x0)²+y²]} + sin[(k0/2(z0+Δz0))(x²+y²)]
% a) 全息图缩放 M（(4.2-3)/(4.2-4)）：H(Mx,My) 中各项变为距离 z/M² 的 FZP、
%    横向位置 x0/M -> M_lat = 1/M, M_long = 1/M²。
% b) 波长缩放 m（(4.2-5)）：重建波长 λ0' = m·λ0（波数 k0/m）-> 焦距 z -> z0/m（m>1 压缩），
%    M_lat = 1, M_long = 1/m。
% c) 无畸变条件：m = 1/M -> M_lat = M_long（Gabor 1949 的原始想法）。
%
% 数值验证（Python 版实跑结论，见笔记图 sim-magnification.png）：
%   a) M = 0.5：缩放全息图在 z0/M² = 120 mm 聚焦出两点，横向间距 = x0/M = 2.00 mm；
%      轴上深度响应聚焦平台中心 120/205 mm（理论 120/200）-> 间距 ≈ Δz0/M² = 80 mm = 4·Δz0 ✓
%   b) m = 2：原始全息图在 z0/m = 50 mm 聚焦，横向位置不变（M_lat = 1）✓
%
% 采样与混叠账（ch3 坑的延续，两组实验用不同 z0 的原因）：
%   FZP 条纹最细周期（视场边缘）= λ0·z0/x_max，需 > 2·dx：
%     a) 缩放全息图等效周期 = λ0·z0/(M²·x_max) -> z0 = 30 mm 即满足；
%     b) 原始全息图 -> 需 z0 >= 2·dx·x_max/λ0 ≈ 81 mm，故取 100 mm。
%   重建传递函数 exp[+jz'k²/2k'] 的相邻格点相位跳变 z'·k·Δk/k' < π：
%     深度扫描全程固定窗（孔径不随 z' 变化），窗缘 12 次方软滚降。
%
% 全书最大的坑：书的前向变换是 +j 指数 -> 前向 ifft2*N²、逆变换 fft2/N²（文件末尾封装）。
% MATLAB 1-based：中心像素 (N/2+1) <-> Python N//2。

clear
lam = 0.6328e-6;  k0 = 2*pi/lam;
x0 = 1.0e-3;                        % 前两点横向间距 (m)
N = 512;  dx = 10e-6;               % 全息图网格
M = 0.5;  m = 2.0;                  % a) 缩放因子  b) 波长缩放因子

x = ((0:N-1) - N/2) * dx;
[X, Y] = meshgrid(x, x);

dkf = 2*pi/(N*dx);                        % 角频率格距 (rad/m)
kx1d = ((0:N-1) - N/2) * dkf;
[KY, KX] = meshgrid(kx1d, kx1d);
KR = hypot(KX, KY);  KR2 = KX.^2 + KY.^2;
k_nyq = 0.55*pi/dx;                       % 通用窗上限

make_h3p = @(z0, dz0) sin(k0*(X.^2+Y.^2)/(2*z0)) + ...    % 三点物体正弦全息图 (4.2-2)
    sin(k0*((X-x0).^2+Y.^2)/(2*z0)) + ...
    sin(k0*(X.^2+Y.^2)/(2*(z0+dz0)));

% ================= a) 全息图缩放 M = 0.5 =================
z0a = 0.030;  dz0a = 0.020;         % 缩放全息图等效条纹周期 λ0z0/(M²x_max) 采样允许更近
Hs = sin(k0*M^2*(X.^2+Y.^2)/(2*z0a)) + ...                % (4.2-4) 逐项解析代入
     sin(k0*M^2*((X-x0/M).^2+Y.^2)/(2*z0a)) + ...
     sin(k0*M^2*(X.^2+Y.^2)/(2*(z0a+dz0a)));
FHs = ft2_book(Hs);

fprintf('--- a) 全息图缩放 M = 0.5（理论：M_lat = 1/M = 2, M_long = 1/M^2 = 4） ---\n');
zs_scan = linspace(0.08, 0.26, 37);
k_win_fix = min(k_nyq, 0.8*pi*k0/(zs_scan(end)*dkf));     % 深度扫描全程同一窗

r1 = recon(FHs, z0a/M^2, k0, k_win_fix, KR, KR2);
[px1, py1] = peaks_of(r1, 0.5, dx);
fprintf('缩放全息图 @ z0/M^2 = %.0f mm：%d 个峰，两点横向间距 = %.2f mm（理论 x0/M = %.2f mm）\n', ...
        z0a/M^2*1e3, numel(px1), px1(end)-px1(1), x0/M*1e3);

I_on = zeros(size(zs_scan));                % 轴上中心像素强度（列最大会被离焦环纹抬高）
for ii = 1:numel(zs_scan)
    rec = recon(FHs, zs_scan(ii), k0, k_win_fix, KR, KR2);
    I_on(ii) = rec(N/2+1, N/2+1);
end
above = I_on > 0.95*max(I_on);              % 聚焦平台（宽度 = 该孔径的景深 ~ ±20 mm）
d = diff([0 above 0]);  st = find(d==1);  en = find(d==-1)-1;
keep = (zs_scan(en)-zs_scan(st)) >= 3*(zs_scan(2)-zs_scan(1));   % 滤掉纹波碎片
st = st(keep);  en = en(keep);
fprintf('轴上深度响应聚焦平台（I>0.95max）：');
for ii = 1:numel(st)
    fprintf('%d~%d mm  ', round(zs_scan(st(ii))*1e3), round(zs_scan(en(ii))*1e3));
end
fprintf('（理论 %d / %d mm）\n', round(z0a/M^2*1e3), round((z0a+dz0a)/M^2*1e3));
cent = 0.5*(zs_scan(st)+zs_scan(en));
fprintf('平台中心间距 = %.0f mm（理论 Δz0/M² = %.0f mm = %.0f·Δz0）\n', ...
        (cent(end)-cent(1))*1e3, dz0a/M^2*1e3, 1/M^2);

% ================= b) 波长缩放 m = 2 =================
z0b = 0.100;  dz0b = 0.030;         % 原始全息图边缘条纹采样要求 z0 >= 81 mm
H3p = make_h3p(z0b, dz0b);
FH3p = ft2_book(H3p);

fprintf('--- b) 波长缩放 m = 2（重建波数 k0/m，理论：M_lat = 1, 焦深 z -> z0/m） ---\n');
r2 = recon(FH3p, z0b/m, k0/m, min(k_nyq, 0.8*pi*(k0/m)/((z0b/m)*dkf)), KR, KR2);
r3 = recon(FH3p, z0b, k0, min(k_nyq, 0.8*pi*k0/(z0b*dkf)), KR, KR2);
px2 = peaks_of(r2, 0.5, dx);
px3 = peaks_of(r3, 0.5, dx);
fprintf('m=2 @ z0/m=%.0f mm 峰位置：%.2f mm；m=1 @ z0=%.0f mm 峰位置：%.2f mm（M_lat = 1）\n', ...
        z0b/m*1e3, px2(2), z0b*1e3, px3(2));
fprintf('第三点焦深：(z0+Δz0)/m = %.0f mm (m=2) vs %.0f mm (m=1) -> 纵向压缩 1/m\n', ...
        (z0b+dz0b)/m*1e3, (z0b+dz0b)*1e3);

% ---- 可视化（对应笔记图 sim-magnification.png） ----
figure(1)
subplot(2,3,1); imagesc(x*1e3, x*1e3, Hs); axis image off; colormap gray
title(sprintf('scaled hologram H(Mx,My), M=%.1f (4.2-3)', M))
subplot(2,3,2); imagesc(x*1e3, x*1e3, r1.^0.4); axis image off
title(sprintf('recon @ z0/M^2=%.0f mm: 2 spots, dx=x0/M', round(z0a/M^2*1e3)))
subplot(2,3,3); plot(zs_scan*1e3, I_on/max(I_on), 'b-', 'LineWidth', 1.2); hold on; grid on
for ii = 1:numel(st)
    xline(zs_scan(st(ii))*1e3, ':b'); xline(zs_scan(en(ii))*1e3, ':b');
end
xline(z0a/M^2*1e3, '--r'); xline((z0a+dz0a)/M^2*1e3, '--r');
xlabel("reconstruction depth z' (mm)"); ylabel('normalized intensity')
title('depth response of scaled hologram (M_{long}=1/M^2)')
subplot(2,3,4); imagesc(x*1e3, x*1e3, H3p); axis image off
title('H_3p(x,y): sine hologram (4.2-2), z0=100 mm')
subplot(2,3,5); imagesc(x*1e3, x*1e3, r2.^0.4); axis image off
title(sprintf('recon @ z0/m=%.0f mm with k''=k0/m (M_{lat}=1)', round(z0b/m*1e3)))
subplot(2,3,6); imagesc(x*1e3, x*1e3, r3.^0.4); axis image off
title(sprintf('recon @ z0=%.0f mm with k0 (reference)', round(z0b*1e3)))

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function F = ft2_book(U)
% 书式前向变换（+j 指数）：F = fftshift(ifft2(ifftshift(U)))*N^2
n = size(U, 1);
F = fftshift(ifft2(ifftshift(U))) * n^2;
end

function U = ift2_book(A)
% 书式逆变换（-j 指数）：U = fftshift(fft2(ifftshift(A)))/N^2
n = size(A, 1);
U = fftshift(fft2(ifftshift(A))) / n^2;
end

function rec = recon(Hspec, zp, k, k_win, KR, KR2)
% 开环重建：|ift2{ Hspec · T(z') }|² 归一化。
% T = 与 h*(z') 对应的频域传递（= conj(OTFosh)，见 ch3 式 3.5-2a）：
% 相邻格点相位跳变 z'·k·Δk/k' < π -> 圆窗半径 k_win；窗缘 12 次方软滚降（硬圆孔有衍射纹波）。
W = 1 ./ (1 + (KR/k_win).^12);
T = W .* exp(1j * zp * KR2 / (2 * k));
rec = abs(ift2_book(Hspec .* T)).^2;
rec = rec / max(rec(:));
end

function px = peaks_of(rec, thr, dx)
% 阈值以上的局部极大峰位置（mm，相对网格中心），按 x 排序；py 未用省略
bw = (imdilate(rec, ones(31)) == rec) & (rec > thr);
CC = bwconncomp(bw);
px = zeros(CC.NumObjects, 1);
for ii = 1:CC.NumObjects
    idx = CC.PixelIdxList{ii};
    [~, im] = max(rec(idx));
    [rr, cc] = ind2sub(size(rec), idx(im));
    px(ii) = (cc - size(rec, 2)/2) * dx * 1e3;
end
px = sort(px);
end
