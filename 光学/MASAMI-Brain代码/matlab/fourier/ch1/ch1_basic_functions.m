% ===== 书 1.5 节：基本函数（复现原书图 1.1 / 图 1.2） =====
% 这些函数是全书建模"孔径 / 光场"的积木，定义必须与原书表 1.2 一致。
% MATLAB/Python 差异：MATLAB 索引从 1 开始（numpy 从 0 开始）；
% MATLAB meshgrid 与 np.meshgrid 默认 'xy' 索引行为一致，圆对称函数下无区别。
% 注：脚本内含局部函数，需要 R2016b 及以上版本。
clear; clc; close all;

% ---------- 自检：定义是否与表 1.2 一致 ----------
fprintf('rect(0)   = %.1f (应为 1.0)\n', rect(0));
fprintf('rect(0.5) = %.1f (应为 0.5)\n', rect(0.5));
fprintf('tri(0.5)  = %.1f (应为 0.5)\n', tri(0.5));
fprintf('sinc(0)   = %.1f (应为 1.0)\n', sincn(0));
fprintf('sinc(0.5) = %.10f (应为 2/pi = %.10f)\n', sincn(0.5), 2/pi);

% ---------- 复现原书图 1.1：一维基本函数 ----------
x  = linspace(-1, 1, 801);
xw = linspace(-4, 4, 1201);
xg = linspace(-2, 2, 801);
figure('Position', [100 100 800 850]);
subplot(3,2,1); plot(x, rect(x));   title('rect(x)');   grid on;
subplot(3,2,2); plot(xw, sincn(xw)); title('sinc(x)');  grid on;
subplot(3,2,3); plot(x, tri(x));    title('\Lambda(x)'); grid on;
subplot(3,2,4);
xn = -2:2;                          % comb 是 delta 周期阵列，画成箭头示意
stem(xn, ones(size(xn)), 'fill');   title('comb(x)'); axis([-2.5 2.5 0 1.2]);
subplot(3,2,5); plot(xg, gauss(xg)); title('exp(-\pix^2)'); grid on;
sgtitle('Basic 1D functions (reproduction of book Fig. 1.1)');
saveas(gcf, 'out_1d_functions.png');

% ---------- 复现原书图 1.2：二维 circ 与高斯 ----------
N  = 80;
xs = ((0:N-1) - N/2) * (4/N);
[X, Y] = meshgrid(xs, xs);          % 'xy' 行为，与 np.meshgrid 默认一致
C  = circle(X, Y, 0.5);             % 表 1.2：半径 1/2
G2 = exp(-pi*(X.^2 + Y.^2));        % 二维高斯（表 1.3 倒数第二行取 a=b=1）
figure;
subplot(1,2,1); mesh(X, Y, C);  title('circ(\surd(x^2+y^2))');
subplot(1,2,2); mesh(X, Y, G2); title('exp[-\pi(x^2+y^2)]');
sgtitle('Examples of 2D functions (reproduction of book Fig. 1.2)');
saveas(gcf, 'out_2d_functions.png');
fprintf('已保存 out_1d_functions.png / out_2d_functions.png\n');

% ---------- 局部函数定义（表 1.2） ----------
function r = rect(t)
    % rect：|t|<1/2 为 1，|t|=1/2 为 1/2，其余 0
    r = zeros(size(t));
    r(abs(t) < 0.5)  = 1;
    r(abs(t) == 0.5) = 0.5;         % 均匀网格很少恰好命中边界
end

function t = tri(u)
    % 三角形函数 Λ：|u|<=1 为 1-|u|，其余 0
    t = zeros(size(u));
    m  = abs(u) <= 1;
    t(m) = 1 - abs(u(m));
end

function s = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)，u=0 处极限 1
    % （与 Signal Processing Toolbox 的 sinc、numpy.sinc 定义一致）
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end

function g = gauss(u)
    % 高斯 exp(-pi*u^2)：Fourier 变换恰好是自己
    g = exp(-pi*u.^2);
end

function c = circle(X, Y, R)
    % 圆盘：sqrt(X^2+Y^2) < R 内取 1（表 1.2 取 R=1/2）
    % 注意：表 1.3 的变换对对应半径 a 的圆盘，与表 1.2 的定义差 2 倍半径，
    % 详见 ch1_verify_circ_jinc.m 的数值验证
    c = double(sqrt(X.^2 + Y.^2) < R);
end
