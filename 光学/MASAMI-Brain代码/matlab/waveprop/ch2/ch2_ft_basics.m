%% ===== 书 Schmidt PM199 第 2 章 =====
%  DFT 基础：定标、重排、正逆变换配对（对应原书 2.1 节与 Listings 2.1-2.4）
%  与 python/waveprop/ch2/ch2_ft_basics.py 数值行为一致（本机无 MATLAB 时
%  以 Python 实跑结果为准，见笔记页面打印输出）。
%  MATLAB 注意：数组从 1 起；偶数 N 时 fftshift 与 ifftshift 相同。
%  需要 R2016b+（脚本内含本地函数）。
clear; clc;

%% 验证 1：式 (2.6) 的 DFT 就是重排后的黎曼和
N = 8; delta = 0.5;
n = (0:N-1) - N/2;                       % n = -N/2 .. N/2-1
x = n * delta;
g = exp(-pi*(x/2).^2);                   % 测试函数（居中排列）
m = (0:N-1) - N/2;
G_sum = zeros(1,N);
for mi = 1:N                             % 逐点按式 (2.4) 的对称指数形式
    G_sum(mi) = delta * sum(g .* exp(-2i*pi*m(mi)*n/N));
end
G_ft = ft(g, delta);
fprintf('验证 1：max|G_ft - G_黎曼和| = %.2e\n', max(abs(G_ft - G_sum)));

%% 验证 2：式 (2.5) 重排公式（原书笔误）
% 原书 (2.5) 第一支 g[(n'+N/2)*delta] 在 n'=1 时坐标 (1+N/2)*delta 越界；
% 第二支 g[(n'-N-2)*delta] 末端只到 -2*delta。正确重排（与 fftshift 一致）：
%   g_tilde{n'} = g[(n'-1)*delta]    n'=1..N/2； g[(n'-1-N)*delta]  n'=N/2+1..N
g_correct = [g(N/2+1:N), g(1:N/2)];      % 即 fftshift(g)（偶数 N）
g_shift = fftshift(g);
fprintf('验证 2：修正公式 vs fftshift 最大差 = %.2e；', max(abs(g_correct - g_shift)));
fprintf('原书第一支越界样本比例 = %.0f%%\n', 100*mean((n + N/2)*delta > x(end)+1e-12));

%% 验证 3：Listing 2.3 高斯 DFT vs 解析 FT（L=5, N=32, a=1）
L = 5; N = 32; a = 1;
delta = L / N;
[n_c, f] = coords(L, N);                 % 居中坐标 x = n*delta 与频率轴
x = n_c * delta;
g_samp = exp(-pi*(a*x).^2);
G_dft = ft(g_samp, delta);
G_ana = 1/a * exp(-pi*(f/a).^2);
[~, i0] = min(abs(f));
fprintf('验证 3：max|G_dft - G_ana| = %.2e（f=0 处 %.2e）\n', ...
    max(abs(G_dft - G_ana)), abs(G_dft(i0) - G_ana(i0)));

%% 验证 4：Listing 2.4 平移高斯的相位（L=5/N=32 与 L=10/N=64 对比）
for pr = [struct('L',5,'N',32), struct('L',10,'N',64)]
    delta = pr.L / pr.N;
    x0 = 5 * delta;
    [n_c, f] = coords(pr.L, pr.N);
    x = n_c * delta;
    G_dft = ft(exp(-pi*(x - x0).^2), delta);
    G_ana = exp(-2i*pi*x0*f) .* exp(-pi*f.^2);
    amp_ok = abs(G_dft) > 1e-3 * max(abs(G_dft));
    ph_err = angle(G_dft(amp_ok) ./ G_ana(amp_ok));
    fprintf('验证 4：L=%-4.1f N=%d：幅度显著样本相位误差 max = %.2e rad\n', ...
        pr.L, pr.N, max(abs(ph_err)));
end

%% 验证 5：ft/ift 往返 + Parseval
L = 5; N = 64;
delta = L / N; delta_f = 1/(N*delta);
[n_c, ~] = coords(L, N);
x = n_c * delta;
g = exp(-pi*(x/1.5).^2) .* (1 + 0.3*cos(2*pi*0.7*x));
g_rt = ift(ft(g, delta), delta_f);
lhs = sum(abs(g).^2) * delta;
rhs = sum(abs(ft(g, delta)).^2) * delta_f;
fprintf('验证 5：max|ift(ft(g))-g| = %.2e；Parseval 相对偏差 %.2e\n', ...
    max(abs(g_rt - g)), abs(lhs-rhs)/lhs);

%% ===== 本地函数 =====
function G = ft(g, delta)
    % 一维 DFT（原书 Listing 2.1）：定标因子 delta 由用户补乘
    G = fftshift(fft(fftshift(g))) * delta;
end

function g = ift(G, delta_f)
    % 一维 DIFT（原书 Listing 2.2）：定标因子 N*delta_f
    g = ifftshift(ifft(ifftshift(G))) * length(G) * delta_f;
end

function [x, f] = coords(L, N)
    % 居中采样坐标（书式 2.21/2.22 同款）：x = n*delta，f = n/(N*delta)，
    % n = -N/2 .. N/2-1。返回升序数组；喂给 ft 前数据须按居中排列。
    delta = L / N;
    n = (-N/2 : N/2-1);
    x = n * delta;
    f = n / (N * delta);
end
