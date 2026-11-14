%% ===== 书 Schmidt PM199 第 1 章（衔接第 2 章主题）=====
%  离散傅里叶变换的三步演化：空间采样 -> 有限截断 -> 频域采样。
%  演示 g(x)=exp(-a|x|) 从连续傅里叶变换到 DFT 的每一步在两个域里各损失/改变了什么。
%  重写自用户旧稿 Alleviating_Effects_of_Discretization.m / 采样结果.py，参数一致：
%  a=10, delta=0.0375, L=0.6, N=16。
clear; clc; close all;

%% 1. 参数
a     = 10;               % 指数衰减系数
delta = 0.0375;           % 空间采样间隔
L     = 0.6;              % 截断窗口宽度
N     = round(L/delta);   % 采样点数 N = 16

%% 2. 连续参考
N_cont = 2000;
x_cont = linspace(-1.5, 1.5, N_cont);
g_cont = exp(-a*abs(x_cont));
f_cont = linspace(-25, 25, N_cont);
G_cont = 2*a ./ (a^2 + (2*pi*f_cont).^2);   % FT{exp(-a|x|)} 解析式

%% 3. 三步演化
% --- 步骤 1：仅采样（空间域乘 comb，频域卷积 comb -> 频谱以 1/delta 为周期复制）---
x_samp = -1.5 : delta : 1.5;
g_samp = exp(-a*abs(x_samp));
G_repl = zeros(size(f_cont));
for p = -2:2                                 % 叠 5 个相邻副本示意混叠
    G_repl = G_repl + 2*a ./ (a^2 + (2*pi*(f_cont - p/delta)).^2);
end

% --- 步骤 2：有限截断（空间域乘 rect(x/L)，频域卷积 L*sinc(L f) -> 波纹）---
x_trunc = (-N/2 : N/2-1) * delta;
g_trunc = exp(-a*abs(x_trunc));
M = 4096;                                    % 高倍补零，模拟"连续的带波纹频谱"
g_pad = zeros(1, M);
g_pad(M/2 - N/2 + 1 : M/2 + N/2) = g_trunc;
G_rip = abs(fftshift(fft(fftshift(g_pad)))) * delta;
f_dense = (-M/2 : M/2-1) / (M*delta);

% --- 步骤 3：频域采样（最终 DFT：频域乘 comb(L f)，空间域周期延拓）---
df = 1/(N*delta);
f_dft = (-N/2 : N/2-1) * df;
G_dft = abs(fftshift(fft(fftshift(g_trunc)))) * delta;

%% 4. 数值验证
fprintf('===== DFT 三步演化数值验证 =====\n');
% (1) 混叠量：步骤 1 的复制谱在 f = 1/(2 delta) 处相对连续谱的抬升
f_nyq = 1/(2*delta);
[~, iny] = min(abs(f_cont - f_nyq));
fprintf('(1) 奈奎斯特频率 %.2f 处 复制谱/连续谱 = %.3f（采样间隙造成的频谱抬升）\n', ...
    f_nyq, G_repl(iny)/G_cont(iny));

% (2) 最终 DFT vs 连续谱：在 DFT 频点 f_m = m/(N delta) 上比较
G_ana_bin = 2*a ./ (a^2 + (2*pi*f_dft).^2);
relerr = max(abs(G_dft - G_ana_bin))/max(G_ana_bin);
fprintf('(2) DFT 谱线 vs 连续解析谱（同频点）最大偏差 = %.2e（截断丢尾 + 邻副本混叠，无法归零）\n', relerr);

% (3) Parseval 能量守恒（DFT 对自身的精确性质，检查实现无误）
lhs = sum(g_trunc.^2) * delta;
rhs = sum(abs(fft(fftshift(g_trunc))*delta).^2) * df;
fprintf('(3) Parseval：sum|g|^2*delta = %.6e，sum|G|^2*df = %.6e，相对偏差 %.2e\n', ...
    lhs, rhs, abs(lhs-rhs)/lhs);

%% 5. 阶梯式对比绘图
figure('Name','DFT 三步演化','Position',[80 50 1200 900]);

subplot(3,2,1);
plot(x_cont, g_cont, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5); hold on;
stem(x_samp, g_samp, 'r', 'Marker', 'none');
xlim([-1.2 1.2]); ylim([0 1.1]); grid on;
title('步骤 1 空间采样：g(x)·comb(x/\delta)（无限采样列）'); ylabel('g(x)');

subplot(3,2,2);
plot(f_cont, G_cont, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5); hold on;
plot(f_cont, G_repl, 'r', 'LineWidth', 2);
xlim([-20 20]); ylim([0 0.25]); grid on;
title('步骤 1 频域：卷积 comb(\deltaf) → 频谱周期复制（混叠）'); ylabel('G(f)');

subplot(3,2,3);
out = abs(x_samp) >= L/2;
stem(x_samp(out), g_samp(out), 'Color', [0.8 0.8 0.8], 'Marker', 'none'); hold on;
stem(x_trunc, g_trunc, 'b', 'Marker', 'o', 'LineWidth', 1.2);
xlim([-1.2 1.2]); ylim([0 1.1]); grid on;
title('步骤 2 有限截断：·rect(x/L)（窗口外采样点丢弃）'); ylabel('g(x)');

subplot(3,2,4);
plot(f_cont, G_repl, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); hold on;
plot(f_dense, G_rip, 'b', 'LineWidth', 2);
xlim([-20 20]); ylim([0 0.25]); grid on;
title('步骤 2 频域：卷积 L·sinc(Lf) → 吉布斯波纹'); ylabel('G(f)');

subplot(3,2,5);
plot(x_cont - L, g_cont, '--', 'Color', [0.75 0.75 0.75]); hold on;
plot(x_cont, g_cont, 'k', 'LineWidth', 1.5);
plot(x_cont + L, g_cont, '--', 'Color', [0.75 0.75 0.75]);
stem(x_trunc, g_trunc, 'g', 'Marker', 'o', 'LineWidth', 1.2);
xlim([-1.2 1.2]); ylim([0 1.1]); grid on;
title('步骤 3 频域采样：空间域卷积 comb(x/L) → 周期延拓'); ylabel('g_t(x)');

subplot(3,2,6);
plot(f_dense, G_rip, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5); hold on;
stem(f_dft, G_dft, 'g', 'Marker', 'o', 'LineWidth', 1.5);
xlim([-20 20]); ylim([0 0.25]); grid on;
title('步骤 3 频域：×comb(Lf) → 最终 N 点 DFT 谱线'); ylabel('G_t(f)');
