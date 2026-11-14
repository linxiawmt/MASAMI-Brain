%% ===== 书 Schmidt PM199 第 2 章 =====
%  二维 DFT：ft2 / ift2 与方形孔径的频谱（对应原书 2.6 节 Listings 2.5-2.6）
%  与 python/waveprop/ch2/ch2_ft2_2d.py 一致（数值以 Python 实跑为准）。
%  采样约定：胞元中心 x=(n+1/2)*delta，孔径边界恰落在胞元边界上，
%  离散孔径面积严格 = D^2（节点采样会把边界样本整格计入）。
%  MATLAB 注意：imagesc + axis xy 才是数学坐标（y 向上）；需要 R2016b+。
clear; clc; close all;

%% 网格与孔径
N = 128; L = 4.0; D = 0.5;
delta = L/N; delta_f = 1/(N*delta);
n = -N/2 : N/2-1;
x = (n + 0.5) * delta;                   % 胞元中心采样
f = n / (N * delta);
[X, Y] = meshgrid(x, x);                 % 与 numpy meshgrid 默认 indexing='xy' 同义
[FX, ~] = meshgrid(f, f);

g = double(abs(X) <= D/2 & abs(Y) <= D/2);
G = ft2(g, delta);

%% 数值验证
fprintf('===== 二维 DFT 验证（方孔 D=%.2f m）=====\n', D);
g_rt = ift2(G, delta_f);
fprintf('(1) max|ift2(ft2(g)) - g| = %.2e\n', max(abs(g_rt(:) - g(:))));
lhs = sum(abs(g(:)).^2) * delta^2;
rhs = sum(abs(G(:)).^2) * delta_f^2;
fprintf('(2) Parseval 相对偏差 = %.2e\n', abs(lhs-rhs)/lhs);
i0 = N/2 + 1;                            % 居中排列的零坐标样本（MATLAB 1-based）
fprintf('(3) |G(0,0)| = %.6f，D^2 = %.6f\n', abs(G(i0,i0)), D^2);
prof = abs(G(i0,:));
prof_ana = abs(D^2 * mysinc(D*f));
ok = abs(f) <= 2/D;
rel = max(abs(prof(ok) - prof_ana(ok))) / D^2;
fprintf('(4) fy=0 剖面 vs 解析 D^2 sinc(D fx)：最大偏差 %.2e（相对峰值）\n', rel);

% 非方阵检查（MATLAB size(G,1) 是行数）：行/列频率间隔不同时，
% 原书 (N*delta_f)^2 定标错误，必须 (M*df_y)*(N*df_x) 分维相乘
M2 = 64;
gx = exp(-pi*x.^2); gy = exp(-pi*x(1:M2).^2);
g2 = gy.' * gx;                          % (M2 x N) 各向异性高斯
G2 = ft2(g2, delta);
rt_book = ifftshift(ifft2(ifftshift(G2))) * (size(G2,1)*delta_f)^2;
rt_fix  = ifftshift(ifft2(ifftshift(G2))) * (M2/(M2*delta)) * (N/(N*delta));
fprintf('(5) %dx%d 非方阵：原书 (行数*df)^2 误差 %.2e；(M*df_y)(N*df_x) 误差 %.2e\n', ...
    M2, N, max(abs(rt_book(:) - g2(:))), max(abs(rt_fix(:) - g2(:))));

%% 绘图
figure('Position',[60 60 1250 400]);
subplot(1,3,1);
imagesc(x, x, g); axis xy equal; colormap gray; colorbar;
title('孔径 g(x,y) = rect(x/D)rect(y/D)'); xlabel('x [m]'); ylabel('y [m]');
subplot(1,3,2);
Gdb = 20*log10(abs(G)/max(abs(G(:)) + eps));
imagesc(f, f, Gdb); axis xy; caxis([-60 0]); colormap parula; colorbar;
title('|G(f_x,f_y)| [dB]，G = D^2 sinc sinc');
xlabel('f_x [m^{-1}]'); ylabel('f_y [m^{-1}]');
subplot(1,3,3);
plot(f(ok), prof_ana(ok), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6); hold on;
plot(f(ok), prof(ok), 'b.', 'MarkerSize', 5);
grid on; xlabel('f_x [m^{-1}]'); ylabel('幅值');
title(sprintf('剖面 f_y=0（最大偏差 %.1e）', rel));
legend({'解析 D^2 sinc(Df_x)', 'DFT 剖面'}, 'FontSize', 9);

%% ===== 本地函数 =====
function G = ft2(g, delta)
    % 二维 DFT（原书 Listing 2.5）
    G = fftshift(fft2(fftshift(g))) * delta^2;
end

function s = mysinc(x)
    % 归一化 sinc：sin(pi x)/(pi x)
    s = ones(size(x));
    nz = x ~= 0;
    s(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end
