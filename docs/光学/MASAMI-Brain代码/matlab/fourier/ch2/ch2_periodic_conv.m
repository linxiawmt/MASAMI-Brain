% ===== 书 2.6 节：周期卷积（式 2.27，复现图 2.7） =====
% 书图 2.7 参数：L = 20, Δx = 1, N = 20 个样本；
%   f = rect，支撑 D1 = 9（样本 x = -4..4 为 1）；h = 右三角形（上升斜坡），
%   支撑 D2 = 15（x ∈ [-7.5, 7.5]，h(x) = (x+7.5)/15）。
% D1 + D2 = 24 > L = 20 → 违反式 (2.27)：FFT 卷积 ≠ 解析卷积（卷绕污染）。
% 演示两条修复路线：(a) 换更大的 L = 40 重新采样；(b) 同一组样本零填充再 FFT。
% MATLAB/Python 差异：conv2/conv 的 'full' 对应 np.convolve 默认；
% 居中坐标写法 n = (0:N-1) - N/2（1-based，零坐标在 N/2+1 号）。
clear; clc; close all;

% ---------- 解析参照：细网格直接卷积 ----------
xa = -30:0.005:30;                          % 中心元素恰为 x=0
dd = xa(2) - xa(1);
conv_full = conv(rectn(xa/9), rampt(xa, 15)) * dd;   % 支撑 [-24, 24]
c0 = (length(conv_full) + 1) / 2;           % conv_full(c0) ↔ t = 0
off = (length(xa) + 1) / 2;                 % xa(off) = 0
conv_ref = conv_full(c0 - off + 1 : c0 + off - 1);   % 对齐到 xa

% ---------- L = 20：违反式 (2.27) ----------
L = 20; N = 20; dx = 1;
n = (0:N-1) - N/2;
x = n * dx;
g = rectn(x/9);                             % D1 = 9
h = rampt(x, 15);                           % D2 = 15
conv_fft = real(fftconv_c(g, h, dx));

err20 = max(abs(conv_ref - conv_fft));
[pkv, ~] = max(conv_ref);
fprintf('L=20（D1+D2=24 > L=20，违反式 2.27）：FFT 卷积与解析最大差异 = %.3f\n', err20);
fprintf('  解析峰值 = %.3f @ x = %.2f；FFT 峰值 = %.3f\n', ...
        conv_ref(pkv), xa(pkv), max(conv_fft));

% ---------- 修复 1：L = 40 重新采样 ----------
L2 = 40; N2 = 40;
n2 = (0:N2-1) - N2/2;
x2 = n2 * dx;
conv_fft2 = real(fftconv_c(rectn(x2/9), rampt(x2, 15), dx));
err40 = max(abs(conv_ref - conv_fft2));
fprintf('L=40（满足式 2.27）：FFT 卷积与解析最大差异 = %.3e\n', err40);

% ---------- 修复 2：零填充（同一组 N=20 样本，填到 N1+N2-1 再 FFT） ----------
P = N + N - 1;                              % 39
g_pad = zeros(1, P); h_pad = zeros(1, P);
o = floor((P - N)/2);
g_pad(o+1:o+N) = g; h_pad(o+1:o+N) = h;
Gp = fft(ifftshift(g_pad)) .* fft(ifftshift(h_pad)) * dx;
conv_pad = real(fftshift(ifft(Gp)));
xp = ((0:P-1) - floor(P/2)) * dx;           % 填充网格的居中坐标
err_pad = max(abs(interp1(xa, conv_ref, xp) - conv_pad));
fprintf('零填充 N1+N2-1=%d：FFT 卷积与解析最大差异 = %.3e\n', P, err_pad);

% ---------- 画图：复现图 2.7(d) + 两种修复 ----------
figure('Position', [60 60 1450 400]);
subplot(1,3,1);
plot(xa, conv_ref, 'k-', 'LineWidth', 1.4); hold on;
plot(x, conv_fft, 'o', 'MarkerSize', 6, 'Color', 'crimson');
grid on; xlim([-30 30]); ylim([-0.3 7]);
xlabel('x'); ylabel('(f*h)(x)');
title(sprintf('L = 20: D1+D2 = 24 > L\nperiodic wrap-around, max err %.2f', err20));
legend('analytic', 'FFT (L=20)', 'FontSize', 8, 'Location', 'northwest');
subplot(1,3,2);
plot(xa, conv_ref, 'k-', 'LineWidth', 1.4); hold on;
plot(x2, conv_fft2, 'o', 'MarkerSize', 5, 'Color', [0 .6 .6]);
grid on; xlim([-30 30]); ylim([-0.3 7]);
xlabel('x'); ylabel('(f*h)(x)');
title(sprintf('L = 40: D1+D2 < L satisfied\nmax err %.1e', err40));
legend('analytic', 'FFT (L=40)', 'FontSize', 8, 'Location', 'northwest');
subplot(1,3,3);
plot(xa, conv_ref, 'k-', 'LineWidth', 1.4); hold on;
plot(xp, conv_pad, 'o', 'MarkerSize', 5, 'Color', [1 .55 0]);
grid on; xlim([-30 30]); ylim([-0.3 7]);
xlabel('x'); ylabel('(f*h)(x)');
title(sprintf('zero-padded to N1+N2-1 = %d\nmax err %.1e', P, err_pad));
legend('analytic', 'zero-padded FFT', 'FontSize', 8, 'Location', 'northwest');
sgtitle('Periodic convolution (book Fig. 2.7): rect D1=9 * right-triangle D2=15');
saveas(gcf, 'out_periodic_conv.png');
fprintf('已保存 out_periodic_conv.png\n');

% ---------- 局部函数 ----------
function r = rectn(t)
    % rect：|t|<1/2 为 1，|t|=1/2 为 1/2，其余 0
    r = zeros(size(t));
    r(abs(t) < 0.5)  = 1;
    r(abs(t) == 0.5) = 0.5;
end

function y = rampt(x, D)
    % 支撑宽度 D 的右三角形（上升斜坡，峰值 1 在右边缘）
    y = zeros(size(x));
    m = abs(x) <= D/2;
    y(m) = (x(m) + D/2) / D;
end

function out = fftconv_c(g, h, dx)
    % FFT 卷积（循环卷积，不填充）：输入居中排列，输出居中排列
    G = fft(ifftshift(g)) .* fft(ifftshift(h)) * dx;
    out = fftshift(ifft(G));
end
