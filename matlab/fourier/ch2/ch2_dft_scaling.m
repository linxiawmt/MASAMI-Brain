% ===== 书 2.3 / 2.5 节：DFT 定标与周期延拓（式 2.19-2.26，复现图 2.6） =====
% 书图 2.6 参数：L = 20, Δx = 1, N = 20 个样本，f(x) = rect(x/5)（半宽 2.5）。
% 验证三件事：
%   (1) 式 (2.19) 的 DFT 定义（核 exp[-j2π(pm/M+qn/N)]）与 MATLAB fft 一致
%       （前提：数据先排成"移位排列"——index 1 对应坐标 0，见 2.4 节）；
%   (2) 乘上 Δx 后，FFT 样本"跟随"周期延拓谱 F_P（Poisson 求和 + 窗平滑），
%       与解析谱 F 的偏差 = 混叠（书 2.5 节：瓣高处样本略大）；
%   (3) 式 (2.19)+(2.20) 连用（1/MN 因子在 ifft 内置）恢复原数组。
% MATLAB/Python 差异：fft/fftshift/ifftshift 语义一致；坐标构造用
% n = (0:N-1) - N/2（1-based 下的居中写法，零坐标在 N/2+1 号）。
clear; clc; close all;
L  = 20.0;        % 采样空间边长 [m]
dx = 1.0;         % 采样间隔 [m]
N  = round(L/dx);
n  = (0:N-1) - N/2;           % 居中整数指标（MATLAB 1-based：坐标 0 在 n=N/2+1 号）
x  = n * dx;                  % x ∈ [-10, 9]，式 (2.21)
f  = n / (N*dx);              % f ∈ [-0.5, 0.4]，式 (2.22)
g  = rectn(x/5);              % rect 半宽 2.5 → 样本 x=-2..2 共 5 个为 1

% ---------- (1) 式 (2.19) 矩阵 DFT vs fft（都在移位排列下） ----------
g_shift = ifftshift(g);                 % 居中 → 移位（index 1 = 坐标 0）
m = 0:N-1;                              % 式 (2.19) 的 m = 0..N-1（移位排列）
p = m.';
W = exp(-2j*pi * p .* m / N);           % DFT 核 exp[-j2π pm/M]
G_matrix = W * g_shift;                 % 式 (2.19)（无 Δx 因子）
G_fft    = fft(g_shift);
err_fft  = max(abs(G_matrix - G_fft));
fprintf('[1] 式(2.19)矩阵 DFT 与 fft 最大差异 = %.2e\n', err_fft);

% ---------- (2) FFT×Δx 样本 = 周期延拓谱 F_P 的样本（近似重合） ----------
G = fftshift(G_fft) * dx;               % 乘 Δx 回连续量纲，再居中排列
F_analytic = 5 .* sincn(5 .* f);        % 解析谱：rect(x/5) <-> 5 sinc(5f)
kk = -3:3;                              % 周期延拓：F_P(f) = Σ_k F(f - k/Δx)
F_periodic = zeros(size(f));
for k = kk
    F_periodic = F_periodic + 5 .* sincn(5 .* (f - k/dx));
end
err_per   = max(abs(G - F_periodic));
err_alias = max(abs(G - F_analytic));
fprintf(['[2] FFT×Δx 与周期谱 F_P 最大差异 = %.2e（空域矩形窗的 sinc 平滑 + ' ...
         '样本化效应，FFT 样本只是"跟随"周期谱）\n'], err_per);
fprintf('    FFT×Δx 与解析谱 F   最大差异 = %.3e（这就是混叠主项）\n', err_alias);
fprintf('    DC 样本 G(0) = %.6f（= 中心纵坐标 ∫g dx = 5）\n', real(G(N/2+1)));

% ---------- (3) 逆变换恢复原数组（式 2.20 的 1/MN 因子在 ifft 内置） ----------
g_back = fftshift(ifft(ifftshift(G / dx)));
err_inv = max(abs(g_back - g));
fprintf('[3] DFT -> DFT^-1 往返恢复原数组最大差异 = %.2e\n', err_inv);

% ---------- 画图（复现书图 2.6） ----------
figure('Position', [80 80 860 720]);
subplot(2,1,1);
plot(x, rectn(x/5), 'k-', 'LineWidth', 1); hold on;
plot(x - L, rectn((x - L)/5), 'k--', 'LineWidth', 0.8);
plot(x + L, rectn((x + L)/5), 'k--', 'LineWidth', 0.8);
plot(x, g, 'o', 'MarkerSize', 6, 'Color', 'crimson');
xlim([-30 30]); ylim([-0.05 1.3]); grid on;
xlabel('x'); ylabel('f(x)');
title('Periodic extension: rect, L = 20, \Deltax = 1');
legend('analytic', 'periodic copies', 'Location', 'northwest');
subplot(2,1,2);
plot(f, abs(F_analytic), 'k-', 'LineWidth', 1.2); hold on;
plot(f, abs(F_periodic), 'k--', 'LineWidth', 0.9);
plot(f, abs(G), 'o', 'MarkerSize', 6, 'Color', 'crimson');
xlim([-1.5 1.5]); ylim([0 5.6]); grid on;
xlabel('f_X'); ylabel('|F(f_X)|');
title('Spectrum: FFT samples follow the periodic spectrum (aliasing lifts the sidelobes)');
legend('analytic F', 'periodic F_P', 'FFT \times \Deltax (N=20)');
saveas(gcf, 'out_periodic_ext.png');
fprintf('已保存 out_periodic_ext.png\n');

% ---------- 局部函数 ----------
function r = rectn(t)
    % rect：|t|<1/2 为 1，|t|=1/2 为 1/2，其余 0（书表 1.2）
    r = zeros(size(t));
    r(abs(t) < 0.5) = 1;
    r(abs(t) == 0.5) = 0.5;
end

function s = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)（与 MATLAB Signal Processing Toolbox 一致）
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end
