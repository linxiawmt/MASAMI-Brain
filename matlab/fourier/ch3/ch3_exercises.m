% ===== 书 3.9 节习题 3.1（三角函数）与 3.3（圆函数）参考解答 =====
% 3.1：底半宽 0.1 m 的三角函数，L = 2、M = 200。书提示的编码：
%       T=1-abs(x); ask=abs(x)<=1; out=T.*mask;
%     其中 ask 与 mask 是同一个变量（原书笔误，两个名字混用，本文件统一
%     写成 mask）。按 3.1(b) 的尺度（底半宽 w = 0.1）缩放后，
%     解析频谱：tri(x/w) <-> w sinc^2(w fx)，w = 0.1 -> 0.1 sinc^2(0.1 fx)。
% 3.3：半径 0.015 m 的圆盘，L = 0.2、M = 200（dx = 0.001，半径恰 15 个样本）。
%     频谱 = jinc 形：R*J1(2*pi*R*rho)/rho，峰值 = 面积 pi*R^2 ≈ 7.07e-4。
%     （半径约定见第 1 章勘误 1：本解按"半径 = R"的 Goodman 约定，
%      与表 1.3 的 circ 变换对一致。）
clear; clc; close all;

% ---------- 习题 3.1：三角函数 ----------
L  = 2; M = 200; dx = L/M;
x  = -L/2:dx:L/2-dx;
w    = 0.1;                     % 底半宽 (m)
T    = 1 - abs(x)/w;            % 书提示第 1 行（按习题尺度缩放）
mask = abs(x) <= w;             % 书提示第 2 行（原书写成 ask——笔误）
tri1 = T.*mask;                 % 书提示第 3 行：边界 |x|=w 处 T=0，无半值歧义

TRI = fftshift(fft(ifftshift(tri1)))*dx;   % 居中 FFT × dx（ifftshift 配对更稳）
fx  = -1/(2*dx) : 1/L : 1/(2*dx) - 1/L;
TRI_an = w*sincn(w*fx).^2;                 % tri(x/w) <-> w sinc^2(w fx)

fprintf('[3.1] 三角函数 FFT 峰值 = %.6f（解析 w = 0.1 = 三角形面积）\n', max(abs(TRI)));
fprintf('[3.1] FFT 与解析最大差异（|fx|<=30 内）= %.2e\n', ...
        max(abs(abs(TRI(abs(fx) <= 30)) - TRI_an(abs(fx) <= 30))));

% ---------- 习题 3.3：圆函数 ----------
L2 = 0.2; M2 = 200; dx2 = L2/M2;
x2 = -L2/2:dx2:L2/2-dx2; y2 = x2;
[X, Y] = meshgrid(x2, y2);
R = 0.015;                       % 圆半径 (m) = 15 个样本
circ2 = circ(sqrt(X.^2+Y.^2)/R); % 调用习题 3.3(a) 建的 circ.m（同文件夹）

C = fftshift(fft2(fftshift(circ2)))*dx2^2; % 二维 FFT × dx^2（偶数 M 下 fftshift 配对成立）
fx2 = -1/(2*dx2) : 1/L2 : 1/(2*dx2) - 1/L2;

CC   = C(M2/2+1, :);                       % fy=0 剖面
J_an = zeros(size(fx2));
nz   = fx2 ~= 0;
J_an(nz)  = R.*besselj(1, 2*pi*R*fx2(nz)) ./ fx2(nz);   % R J1(2 pi R rho)/rho
J_an(~nz) = pi*R^2;                                     % rho=0 极限 = 面积

fprintf('[3.3] 圆盘 FFT 峰值 = %.6e（解析面积 pi*R^2 = %.6e）\n', ...
        max(abs(C(:))), pi*R^2);
fprintf('[3.3] fy=0 剖面与 jinc 最大差异（|fx|<=30 内）= %.2e\n', ...
        max(abs(abs(CC(abs(fx2) <= 30)) - abs(J_an(abs(fx2) <= 30)))));

% ---------- 画图 ----------
figure('Position', [60 60 1150 700]);
subplot(2,2,1);
plot(x, tri1, '-o', 'MarkerSize', 4); grid on; axis([-0.3 0.3 0 1.2]);
xlabel('x (m)'); title('Ex 3.1: triangle, half-base 0.1 m');
subplot(2,2,2);
plot(fx, abs(TRI), fx, TRI_an, ':'); grid on; xlim([-60 60]);
xlabel('fx (cyc/m)'); title('Ex 3.1: |FFT| vs w sinc^2(w f_x)');
legend('discrete', 'analytic', 'FontSize', 8);
subplot(2,2,3);
imagesc(x2, y2, circ2); colormap('gray'); axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)'); title('Ex 3.3: circle, R = 0.015 m');
subplot(2,2,4);
plot(fx2, abs(CC), 'o', 'MarkerSize', 4, fx2, abs(J_an), ':'); grid on;
xlim([-60 60]); xlabel('fx (cyc/m)');
title('Ex 3.3: slice vs jinc R J_1(2\pi R\rho)/\rho');
legend('discrete', 'analytic', 'FontSize', 8);
saveas(gcf, 'out_exercises.png');
fprintf('已保存 out_exercises.png\n');

% ---------- 局部函数 ----------
function s = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)，u=0 处极限 1（表 1.2；np.sinc 同定义）
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end
