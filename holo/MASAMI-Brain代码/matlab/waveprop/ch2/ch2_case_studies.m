%% ===== 书 Schmidt PM199 第 2 章 =====
%  三个变换案例研究：sinc / 高斯 / 带二次相位的高斯（对应原书 2.5 节）
%  与 python/waveprop/ch2/ch2_case_studies.py 一致（数值以 Python 实跑为准）。
%  含原书勘误复现：2.5.3 节 "f_{x,0.01}=1.96 m^-1" 实为 p=0.001 的值，
%  式 (2.31) 代入 p=0.01 应为 1.60 m^-1。需要 R2016b+。
clear; clc;

%% 案例 1：sinc 信号（原书图 2.10：a=1.1, delta=0.85/a, N=32）
a = 1.1; delta = 0.85/a; N = 32;
[x, f] = coords(N*delta, N);
g = mysinc(a*x);
G = ft(g, delta);
G_a = ft_ana_sinc(f, a);
band = abs(f) <= a/2;
rip = max(abs(abs(G(band)) - abs(G_a(band)))) * a;    % 相对峰值 1/a
ph_edge = angle(G(~band));
frac_pi = mean(abs(abs(ph_edge) - pi) < 0.1);
fprintf('案例 1：带内幅度波纹 = %.2e；带外样本相位≈pi 比例 = %.0f%%\n', ...
    rip, 100*frac_pi);

%% 案例 2：高斯信号（原书图 2.11/2.12：a=1）
a = 1.0;
f_e2 = a * sqrt(2/pi);
delta_e2 = 1/(2*f_e2);                    % 1/e^2 判据（原书式 2.25）
delta_01 = 1/(2*sqrt(a^2/pi*log(100)));   % p=0.01 判据（原书式 2.27）
for case2 = [struct('tag','delta_e2  ','d',delta_e2), ...
             struct('tag','delta_0.01','d',delta_01)]
    N = 32;
    [x, f] = coords(N*case2.d, N);
    G = ft(exp(-pi*(a*x).^2), case2.d);
    G_left = real(G(1));                  % 居中排列最左样本 = 频带负边缘
    G_anL = 1/a * exp(-pi*f(1)^2/a^2);
    fprintf('案例 2：%s = %.4f：最左样本 %.4f（解析 %.4f，混叠抬升占峰值 %.0f%%）\n', ...
        case2.tag, case2.d, G_left, G_anL, 100*(G_left - G_anL)/(1/a));
end

%% 案例 3：二次相位高斯（原书图 2.13-2.15），a=0.25、b=0.57
a = 0.25; b = 0.57;
for p = [0.01, 0.001]
    fp = bandwidth_p(a, b, p);            % 原书式 (2.31)
    fprintf('案例 3：p=%-6g：式(2.31) f_p = %.4f m^-1%s\n', p, fp, ...
        ternary(p==0.001, '  <- 原书正文的 1.96 实为这一档（勘误）', ''));
end
% 原书采样参数下的 DFT 精度
cfgs = {struct('tag','原书 delta=0.25 m（=p=0.001 判据）','d',0.25,'N',40), ...
        struct('tag','p=0.01 判据 delta=0.312 m            ','d', ...
               1/(2*bandwidth_p(a,b,0.01)),'N',32)};
for c = 1:2
    d = cfgs{c}.d; N = cfgs{c}.N;
    [x, f] = coords(N*d, N);
    g = exp(-pi*(a*x).^2) .* exp(1i*pi*(b*x).^2);
    G = ft(g, d);
    G_a = 1/sqrt(a^2 - 1i*b^2) * exp(-pi*f.^2/(a^2 - 1i*b^2));
    ok = abs(G_a) > 1e-3 * max(abs(G_a));
    amp_err = max(abs(abs(G(ok)) - abs(G_a(ok)))) / max(abs(G_a));
    ph_err = max(abs(angle(G(ok) ./ G_a(ok))));
    fprintf('案例 3：%s N=%d：带内幅度误差 %.2e，相位误差 %.2e rad\n', ...
        cfgs{c}.tag, N, amp_err, ph_err);
end
% 原书图 2.14 的带宽（a=0.33，p=0.01）
fprintf('案例 3：图 2.14 带宽（a=0.33，p=0.01）：');
for bb = [0.75, 1.5, 2.5]
    fprintf('b=%.2f -> %.2f  ', bb, bandwidth_p(0.33, bb, 0.01));
end
fprintf('m^-1\n');

%% ===== 本地函数 =====
function G = ft(g, delta)
    % 一维 DFT（原书 Listing 2.1）
    G = fftshift(fft(fftshift(g))) * delta;
end

function [x, f] = coords(L, N)
    % 居中采样坐标：x = n*delta，f = n/(N*delta)，n = -N/2 .. N/2-1
    delta = L / N;
    n = (-N/2 : N/2-1);
    x = n * delta;
    f = n / (N * delta);
end

function s = mysinc(x)
    % 归一化 sinc：sin(pi x)/(pi x)（避免依赖 Signal Processing Toolbox）
    s = ones(size(x));
    nz = x ~= 0;
    s(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end

function G = ft_ana_sinc(f, a)
    % sinc(a x) 的解析 FT：(1/a) rect(f/a)
    G = (abs(f) < a/2) / a + (abs(f) == a/2) * (0.5/a);
end

function fp = bandwidth_p(a, b, p)
    % p 电平带宽（原书式 2.31；b=0 退化为式 2.27）
    fp = sqrt(-(a^2 + b^4/a^2)/pi * log(p));
end

function out = ternary(cond, vTrue, vFalse)
    if cond, out = vTrue; else, out = vFalse; end
end
