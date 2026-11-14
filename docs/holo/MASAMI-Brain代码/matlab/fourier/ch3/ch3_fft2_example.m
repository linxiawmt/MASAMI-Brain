% ===== 书 3.7 节：二维矩形函数、meshgrid、imagesc 与二维 FFT =====
% 参数：wx = 0.1 m、wy = 0.05 m 的矩形（0.2 × 0.1 m），L = 2 m、M = 200/边
% -> dx = 0.01 m，200×200 阵列对应 2 × 2 m。
% MATLAB/Python 差异：
%   - meshgrid 两边行为一致（'xy'：X 的行是 x 的副本、Y 的列是 y 的副本；
%     Python 若写 indexing='ij' 两维会互换）；
%   - (M/2+1,:) 选第 M/2+1 行 = y=0 的行（行 <-> y！），Python 是 g[M//2, :]；
%   - imagesc 默认 y 向下，axis xy 翻转；Python imshow 用 origin='lower'；
%   - 二维必须 fft2/ifft2（Python np.fft.fft2）：一维 fft 作用于二维阵列
%     只变换行向量，得到"条纹"状错误结果（书 3.8 节原话）。
clear; clc; close all;

wx = 0.1;                   % 矩形 x 半宽 (m)
wy = 0.05;                  % 矩形 y 半宽 (m)
L  = 2;                     % 边长 x&y (m)
M  = 200;                   % 每边样本数
dx = L/M;                   % 采样间隔 (m)
x  = -L/2:dx:L/2-dx;        % x 坐标
y  = x;                     % y 坐标
[X, Y] = meshgrid(x, y);    % X 的行 = x 的副本；Y 的列 = y 的副本
g = rect(X/(2*wx)).*rect(Y/(2*wy));   % 二维矩形（可分离乘积，调 rect.m）

% ---------- 图像显示（书图 3.10a）与中心切片（书图 3.10b） ----------
% imagesc 把数组缩放到全显示范围；colormap('gray') 灰度；axis square 方形
% 边框；axis xy 让 y 自下而上——常规图像文件第 1 行在顶部，MATLAB 图像
% 显示默认 y 向下增，axis xy 把它翻回笛卡尔习惯。
figure('Position', [60 60 1000 420]);
subplot(1,2,1);
imagesc(x, y, g); colormap('gray'); axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)'); title('g(x,y): 2D rect (Fig. 3.10a)');
subplot(1,2,2);
plot(x, g(M/2+1,:), '-o', 'MarkerSize', 4); grid on;   % 第 M/2+1 行 = y=0
axis([-1 1 0 1.2]); xlabel('x (m)');
title('x-slice through g(M/2+1,:) (Fig. 3.10b)');

% ---------- 二维 FFT ----------
g0 = fftshift(g);           % 预移位（偶数 M：fftshift = ifftshift）
G0 = fft2(g0) * dx^2;       % 二维 FFT 并乘 dx^2 定标（式 2.19 的二维版）
G  = fftshift(G0);          % 居中
fx = -1/(2*dx) : 1/L : 1/(2*dx) - 1/L;   % 频率坐标（式 2.22）
fy = fx;

% ---------- 变换幅度显示：曲面 + 中心剖面（书图 3.11） ----------
% surf 曲面显示；camlight/lighting/shading 只是化妆。fy=0 剖面有两条参照：
%   (1) 连续矩形解析谱 4*wx*wy*|sinc(2*wx*fx)|，峰值 = 面积 0.020；
%   (2) 离散精确谱（Dirichlet 核）：|G(fx,0)| = dx^2*ny*|sin(nx*pi*fx*dx)|
%       ./|sin(pi*fx*dx)|，nx/ny 是 x/y 向值为 1 的样本数（21、11）——
%       FFT 剖面应与它逐点重合，峰值 = 离散面积 21*11*dx^2 = 0.0231。
%       书图 3.11(b) 的峰 ~0.023 正是这个离散面积：wx=0.1、wy=0.05 都是
%       "偶数样本宽"，书 rect 的奇数样本约定各多给一圈边界样本（3.1 节），
%       不是混叠。
nx = sum(rect(X/(2*wx))(:));          % x 向值为 1 的样本数（21）
ny = sum(rect(Y/(2*wy))(:));          % y 向值为 1 的样本数（11）
D  = abs(sin(nx*pi*fx*dx) ./ sin(pi*fx*dx));
D(fx == 0) = nx;                      % 0/0 极限
slice_exact = dx^2 * ny * D;
figure('Position', [80 80 1150 420]);
subplot(1,2,1);
surf(fx, fy, abs(G), 'EdgeColor', 'none'); colormap('gray');
xlabel('fx (cyc/m)'); ylabel('fy (cyc/m)');
title('|G(f_x,f_y)| surface (Fig. 3.11a)');
subplot(1,2,2);
plot(fx, abs(G(M/2+1,:)), fx, slice_exact, ':', ...
     fx, 4*wx*wy*abs(sincn(2*wx*fx)), '--'); grid on;
xlabel('fx (cyc/m)'); xlim([-50 50]);
title('fx slice: FFT vs two references (Fig. 3.11b)');
legend('FFT', 'discrete exact (Dirichlet)', '4w_xw_y sinc(2w_x f_x)', 'FontSize', 8);

fprintf('x/y 向值为 1 的样本数 = %d x %d（边界取 1 的奇数样本约定）\n', nx, ny);
fprintf('DC 样本 |G(0,0)| = %.6f = 离散面积 %d*%d*dx^2 = %.4f\n', ...
        abs(G(M/2+1, M/2+1)), nx, ny, nx*ny*dx^2);
fprintf('  （连续矩形面积 4*wx*wy = %.3f；书图 3.11(b) 的峰 ~0.023 就是离散面积）\n', ...
        4*wx*wy);
in_main = abs(fx) <= 25;
fprintf('fy=0 剖面 vs 离散精确谱（Dirichlet）最大差异 = %.2e（机器精度）\n', ...
        max(abs(abs(G(M/2+1,in_main)) - slice_exact(in_main))));
fprintf('fy=0 剖面 vs 连续解析谱最大差异 = %.4f（DC 处 0.0231 vs 0.020 + 零点错位）\n', ...
        max(abs(abs(G(M/2+1,in_main)) - 4*wx*wy*abs(sincn(2*wx*fx(in_main))))));
saveas(gcf, 'out_fft2_example.png');
fprintf('已保存 out_fft2_example.png\n');

% ---------- 局部函数 ----------
function s = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)，u=0 处极限 1（表 1.2；np.sinc 同定义）
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end
