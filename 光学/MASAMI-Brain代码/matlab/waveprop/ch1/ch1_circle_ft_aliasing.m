%% ===== 书 Schmidt PM199 第 1 章 =====
%  圆孔的二维傅里叶变换 = 夫琅禾费衍射（艾里斑），以及粗采样带来的频谱混叠。
%  重写自用户旧稿 Fourier_Transfer_2D.m / circle_digital.m / 二维傅里叶变换.py / 圆孔离散化.py：
%    - 细/粗两套网格合并到一个脚本，补上与连续理论的定量核对；
%    - 统一"胞元中心"采样约定。
%  MATLAB 注意：fftshift 对每维同样适用；imagesc 默认第 1 行在顶部，axis xy 翻转。
clear; clc; close all;

%% 1. 参数
D   = 100e-6;         % 圆孔直径 100 um（circ(2r/D)，Appendix A）
wvl = 532e-9;

% --- A. 细网格：逼近连续圆孔 ---
N_f = 512;
dx_f = 1e-6;                          % 1 um 像素，圆孔直径占 100 像素
x_f = ((1:N_f) - N_f/2 - 0.5) * dx_f;
[X_f, Y_f] = meshgrid(x_f, x_f);
U_f = double(sqrt(X_f.^2 + Y_f.^2) <= D/2);       % circ(2r/D)

% --- B. 粗网格：同一个圆孔，像素放大 4 倍 ---
N_c = 128;
dx_c = 4e-6;                          % 4 um 像素，圆孔直径只占 25 像素 -> 边缘锯齿
x_c = ((1:N_c) - N_c/2 - 0.5) * dx_c;
[X_c, Y_c] = meshgrid(x_c, x_c);
U_c = double(sqrt(X_c.^2 + Y_c.^2) <= D/2);

%% 2. 居中二维傅里叶变换（ft2，含面积元 dx^2）
G_f = ft2(U_f, dx_f);   I_f = abs(G_f).^2;   I_f = I_f/max(I_f(:));
G_c = ft2(U_c, dx_c);   I_c = abs(G_c).^2;   I_c = I_c/max(I_c(:));
f_f = ((1:N_f) - N_f/2 - 1) / (N_f*dx_f);         % 频率轴 (1/m)
f_c = ((1:N_c) - N_c/2 - 1) / (N_c*dx_c);

%% 3. 数值验证
fprintf('===== 圆孔傅里叶变换数值验证 =====\n');
% (1) 直流值：FT{circ(2r/D)}(0) = 圆孔面积 = pi*D^2/4
dc_th  = pi*D^2/4;
dc_num = abs(G_f(N_f/2+1, N_f/2+1));
fprintf('(1) 直流分量理论 pi*D^2/4 = %.4e，数值 = %.4e，相对偏差 %.2e\n', ...
    dc_th, dc_num, abs(dc_num-dc_th)/dc_th);

% (2) 艾里斑第一暗环：J1 零点 3.8317 -> f = 3.8317/(pi*D)；三点抛物线细化极小位置
f1_th = 3.8317/(pi*D);
row  = I_f(N_f/2+1, :);                       % 沿 fx 轴（fy=0）的剖面
half = N_f/2 + 2;                             % f>0 半轴起点（f(0) 在第 N/2+1 个）
idx = find(row(half:end-1) < row(half-1:end-2) & row(half:end-1) < row(half+1:end)) + half - 1;
i0 = idx(1);
df_bin = f_f(2) - f_f(1);
den = row(i0-1) - 2*row(i0) + row(i0+1);
f1_num = f_f(i0) + 0.5*df_bin*(row(i0-1) - row(i0+1))/(den + 1e-30);
fprintf('(2) 第一暗环理论 f = %.3f mm^-1，数值（细化）= %.3f mm^-1，相对偏差 %.2e\n', ...
    f1_th*1e3, f1_num*1e3, abs(f1_num-f1_th)/f1_th);

% (3) 粗网格的两类失真（沿 fx 剖面对比细网格）：
%     像素化边缘把光散射进艾里暗环（低频段暗环被"填亮"）；
%     频谱以 1/dx_c 为周期复制，在近奈奎斯特端（|f|>100 mm^-1）出现混叠尾巴。
fx_c_mm = f_c*1e-3;
I_f_interp = interp1(f_f*1e-3, row/row(half-1), fx_c_mm, 'linear', NaN);
prof_c = I_c(N_c/2+1, :)/I_c(N_c/2+1, N_c/2+1);
[~, i_ring] = min(abs(fx_c_mm - f1_th*1e-3));
fill = prof_c(i_ring)/max(I_f_interp(i_ring), 1e-12);
high = abs(fx_c_mm) > 100;
err_high = max(abs(prof_c(high) - I_f_interp(high)), [], 'omitnan');
fprintf('(3) 第一暗环处 粗/细强度比 = %.1f（像素化把暗环填亮）；|f|>100 mm^-1 复制尾巴最大偏差 = %.2e\n', ...
    fill, err_high);

%% 4. 绘图
figure('Name','圆孔频谱与采样混叠','Position',[60 60 1250 800]);

subplot(2,3,1);
imagesc(x_f*1e6, x_f*1e6, U_f); axis xy; axis image; colormap gray;
title('(a) 细网格圆孔（连续近似）'); xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2,3,2);
imagesc(f_f*1e-3, f_f*1e-3, log10(I_f + 1e-6)); axis xy; axis image;
caxis([-5 0]); colormap gray; colorbar;
title('(b) 艾里斑频谱（对数尺度）'); xlabel('f_x (mm^{-1})'); ylabel('f_y (mm^{-1})');

subplot(2,3,3);
s = 4;                                        % 3D 面抽稀，仅中心 1/4 区域
c = (N_f/2 - 32*s) : s : (N_f/2 + 32*s);
surf(f_f(c)*1e-3, f_f(c)*1e-3, I_f(c, c), 'EdgeColor', 'none'); shading interp;
colormap jet;
title('(c) 艾里斑 3D 光强'); xlabel('f_x (mm^{-1})'); ylabel('f_y (mm^{-1})');
view(-37.5, 30);

subplot(2,3,4);
imagesc(x_c*1e6, x_c*1e6, U_c); axis xy; axis image; colormap gray;
title('(d) 粗网格圆孔（边缘锯齿）'); xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2,3,5);
imagesc(f_c*1e-3, f_c*1e-3, log10(I_c + 1e-6)); axis xy; axis image;
caxis([-5 0]); colormap gray; colorbar;
title('(e) 粗网格频谱：高频端出现混叠'); xlabel('f_x (mm^{-1})'); ylabel('f_y (mm^{-1})');

subplot(2,3,6);
semilogy(f_f*1e-3, row/row(half-1), 'b-', 'LineWidth', 1.2); hold on;
semilogy(fx_c_mm, prof_c, 'r.-', 'MarkerSize', 4);
xline(f1_th*1e-3, ':g', 'LineWidth', 1.2);
xlim([-60 60]); ylim([1e-5 1.1]); grid on;
xlabel('f_x (mm^{-1})'); ylabel('归一化光强');
legend('细网格 (dx=1 \mum)', '粗网格 (dx=4 \mum)', '理论第一暗环 12.2 mm^{-1}', 'Location', 'south');
title('沿 f_x 剖面对比（对数纵轴）');

%% ===== 本地函数 =====
function G = ft2(g, delta)
    % G = fftshift(fft2(fftshift(g)))*delta^2：中心化二维傅里叶变换（连续 FT 的离散近似）
    G = fftshift(fft2(fftshift(g))) * delta^2;
end
