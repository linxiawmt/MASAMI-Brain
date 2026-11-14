% ===== 书 1.5 节：用 FFT 数值验证表 1.3 的一维变换对 =====
% 连续傅里叶变换的 DFT 近似（矩形求和）：
%     G(f) ≈ dx * Σ_n g(x_n) exp(-j2*pi*f*x_n)
% MATLAB fft 与 numpy.fft 一样不对正变换归一化（归一化在逆变换里），
% 因此乘 dx 就回到连续 FT 的量纲。
% 信号定义在 x ∈ [-L/2, L/2)：ifftshift -> fft -> fftshift 得"居中"频谱。
% （偶数 N 时 ifftshift 与 fftshift 相同，这里仍写全流程以对应 Python 代码。）
% 定标与居中的系统讲解见本书第 2、3 章。
clear; clc; close all;
N  = 512;             % 采样点数（偶数）
L  = 8.0;             % 空间窗口宽度 [m]
dx = L / N;           % 采样间隔 [m]
n  = (0:N-1) - N/2;   % MATLAB 1-based 索引下构造居中坐标的常用写法
x  = n * dx;          % x ∈ [-L/2, L/2)
f  = n / (N * dx);    % 频率轴 [cycles/m]，与 x 同序（居中）

% ---------- 对 1：rect(x) <-> sinc(f)（表 1.3 第 1 行，a=1） ----------
G1 = ft_centered(rect(x), dx);
e1 = max(abs(G1 - sincn(f)));
fprintf('[1] rect(x) <-> sinc(f)           最大绝对误差 = %.3e\n', e1);

% ---------- 对 2：Λ(x) <-> sinc^2(f)（表 1.3 第 3 行，a=1） ----------
G2 = ft_centered(tri(x), dx);
e2 = max(abs(G2 - sincn(f).^2));
fprintf('[2] tri(x) <-> sinc^2(f)          最大绝对误差 = %.3e\n', e2);

% ---------- 对 3：exp(-pi x^2) <-> exp(-pi f^2)（表 1.3 第 5 行，a=1） ----------
g3 = exp(-pi*x.^2);
G3 = ft_centered(g3, dx);
e3 = max(abs(G3 - exp(-pi*f.^2)));
fprintf('[3] exp(-pi x^2) <-> exp(-pi f^2)  最大绝对误差 = %.3e\n', e3);

% ---------- 定理 1：相似性 rect(x/2) <-> 2 sinc(2f)（表 1.1，a=2） ----------
G4 = ft_centered(rect(x/2), dx);
e4 = max(abs(G4 - 2*sincn(2*f)));
fprintf('[4] 相似性 rect(x/2) <-> 2sinc(2f) 最大绝对误差 = %.3e\n', e4);

% ---------- 定理 2：位移 rect(x-0.7) <-> sinc(f) exp(-j2*pi*f*0.7) ----------
x0 = 0.7;
G5 = ft_centered(rect(x - x0), dx);
e5 = max(abs(G5 - sincn(f) .* exp(-2j*pi*f*x0)));
fprintf('[5] 位移 rect(x-0.7) 复数最大误差 = %.3e\n', e5);

% ---------- 定理 3：中心纵坐标 G(0) = ∫ g dx（表 1.1 最后一行） ----------
fprintf(['[6] 中心纵坐标 G(0): 数值 = %.6f, 求和*dx = %.6f, 解析 = 1.0\n'], ...
        real(G1(N/2+1)), sum(rect(x))*dx);

% ---------- 勘误验证：表 1.3 最后一行 chirp 变换对（符号有误） ----------
% 原书：exp(-j*pi*x^2) <-> j*exp(-j*pi*f^2)（2D 版本 j*|ab|*exp[-j*pi(...)]）
% 正确：exp(-j*pi*x^2) <-> exp(-j*pi/4)*exp(+j*pi*f^2)
%       （2D 版本 -j*|ab|*exp[+j*pi(a^2 fX^2 + b^2 fY^2)]）
% chirp 不衰减，截断积分误差 ~ 1/(窗口半宽)，用宽窗口就能清晰区分两种公式。
Nc = 4096; Lc = 40.0; dxc = Lc/Nc;
nc = (0:Nc-1) - Nc/2;
xc = nc*dxc; fc = nc/(Nc*dxc);
gc = exp(-1j*pi*xc.^2);
Gc = fftshift(fft(ifftshift(gc)))*dxc;
m  = abs(fc) <= 1;
mine = exp(-1j*pi/4)*exp(1j*pi*fc.^2);
book = 1j*exp(-1j*pi*fc.^2);
fprintf(['[7] chirp 修正式误差 = %.3e, 原书式误差 = %.3e（应前者 << 后者）\n'], ...
        max(abs(Gc(m) - mine(m))), max(abs(Gc(m) - book(m))));

% ---------- 画对比图：数值（圆点） vs 解析（实线） ----------
pairs = {G1, sincn(f), 'rect <-> sinc'; ...
         G2, sincn(f).^2, 'tri <-> sinc^2'; ...
         G3, exp(-pi*f.^2), 'gauss <-> gauss'};
figure('Position', [100 100 1200 380]);
for k = 1:3
    subplot(1, 3, k);
    A = pairs{k, 2}; Gr = real(pairs{k, 1});
    m = abs(f) <= 4;
    plot(f(m), A(m), 'b-', 'LineWidth', 1.5); hold on;
    plot(f(m), Gr(m), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'none');
    grid on; xlabel('f [cycles/m]'); title(pairs{k, 3});
    legend('analytic', 'DFT numeric');
end
sgtitle(sprintf('Table 1.3 pairs verified with FFT (N=%d, L=%.0f)', N, L));
saveas(gcf, 'out_1d_pairs.png');
fprintf('已保存 out_1d_pairs.png\n');

% ---------- 局部函数 ----------
function G = ft_centered(g, dx)
    % 一维连续 FT 的 DFT 近似（居中输出）；dx 需显式传入
    % （MATLAB 局部函数不像 Python 函数那样能闭包捕获脚本变量）
    G = fftshift(fft(ifftshift(g))) * dx;
end

function r = rect(t)
    r = zeros(size(t));
    r(abs(t) < 0.5)  = 1;
    r(abs(t) == 0.5) = 0.5;
end

function t = tri(u)
    t = zeros(size(u));
    m = abs(u) <= 1;
    t(m) = 1 - abs(u(m));
end

function s = sincn(u)
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end
