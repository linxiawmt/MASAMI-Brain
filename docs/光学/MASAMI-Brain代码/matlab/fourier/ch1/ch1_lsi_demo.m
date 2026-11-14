% ===== 书 1.6 节：线性空间不变系统 —— 卷积定理数值演示（式 1.11-1.13） =====
% 线性空间不变系统的输出是卷积：  g2(x,y) = (g1 * h)(x,y)
% 两条等价路线：
%   A. 直接卷积 conv2（对应式 1.11 的积分求和，物理直观）
%   B. 频域法（式 1.13）：FFT(g1) .* H -> IFFT，H 是传递函数
% 本脚本用"零填充 + FFT"实现线性卷积，与 conv2 的直接结果对比，
% 两者应只差浮点舍入（~1e-15）。
% 注意：不填充直接 fft2 相乘得到的是"循环卷积"，边界会互相缠绕——
% 这是第 2 章"计算以外的世界"和第 5 章角谱法的核心话题，这里先埋个伏笔。
% 显示差异：Python 的 imshow 用 origin='lower' 让 y 轴向上；
% MATLAB 的 imagesc 默认 y 轴向下，需 axis xy 达到同样的显示方向。
clear; clc; close all;
N  = 256;           % 物面网格每边点数
L  = 1.0;           % 物面视场 [m]
dx = L / N;
x  = ((0:N-1) - N/2) * dx;
[X, Y] = meshgrid(x, x);

% ---------- 输入物 g1：圆盘 + 小方块 + 双缝（简单"物体"） ----------
g1 = double(sqrt((X+0.25).^2 + (Y+0.22).^2) < 0.10);            % 圆盘
g1 = g1 + double((abs(X-0.22) < 0.05) & (abs(Y+0.20) < 0.05));  % 方块
g1 = g1 + double((abs(X-0.05) < 0.012) & (abs(Y-0.25) < 0.10)); % 细棒 1
g1 = g1 + double((abs(X+0.05) < 0.012) & (abs(Y-0.25) < 0.10)); % 细棒 2
g1 = min(g1, 1);

% ---------- 冲激响应 h：小圆盘（"模糊斑"，即光学里的点扩散函数 PSF） ----------
rh = 0.03;
Kh = 31;                              % h 网格边长（奇数，中心对齐）
xh = ((0:Kh-1) - floor(Kh/2)) * dx;
[Xh, Yh] = meshgrid(xh, xh);
h = double(sqrt(Xh.^2 + Yh.^2) < rh);
h = h / sum(h(:));                    % 归一化能量：直流增益 H(0,0)=1，像不整体变亮

% ---------- 路线 B：频域法（零填充到 N+Kh-1，做线性卷积） ----------
P = N + Kh - 1;
G2full = real(ifft2(fft2(g1, P, P) .* fft2(h, P, P)));
o = floor(Kh/2);
g2_freq = G2full(o+1 : o+N, o+1 : o+N);
% 裁切说明：h 中心在 (0,0)，其样本在填充数组中平移了 floor(Kh/2)；
% 线性卷积结果里 g1 的第 u 个样本对应整体输出的第 u + floor(Kh/2) 个。

% ---------- 路线 A：直接二维卷积（conv2 的 'same' 与 scipy 对齐） ----------
g2_direct = conv2(g1, h, 'same');

diffval = max(abs(g2_direct(:) - g2_freq(:)));
rel = abs(g2_direct - g2_freq) / max(g2_freq(:));   % 相对误差（显示用）
fprintf('两种路线结果最大差异 = %.3e（相对量 %.1e，浮点舍入量级）\n', ...
        diffval, max(rel(:)));

% ---------- 画图 ----------
h_full = zeros(N);                   % 把 h 镶进全画面，统一显示比例
oh = floor((N-Kh)/2);
h_full(oh+1:oh+Kh, oh+1:oh+Kh) = h;
figure('Position', [80 80 1500 380]);
subplot(1,4,1); imagesc(x, x, g1);      axis xy off; title('input g1 (object)');
colormap(gca, gray);
subplot(1,4,2); imagesc(x, x, h_full);  axis xy off;
title(sprintf('impulse response h (PSF), %dx%d grid', Kh, Kh));
colormap(gca, gray);
subplot(1,4,3); imagesc(x, x, g2_freq); axis xy off; title('output g2 = g1 * h');
colormap(gca, gray);
subplot(1,4,4); imagesc(x, x, rel);     axis xy off; colorbar;
title(sprintf('|direct - spectral| / max(g2)\n(max %.0e)', max(rel(:))));
sgtitle('Linear space-invariant system: convolution vs. transfer-function route');
saveas(gcf, 'out_lsi_demo.png');
fprintf('已保存 out_lsi_demo.png\n');
