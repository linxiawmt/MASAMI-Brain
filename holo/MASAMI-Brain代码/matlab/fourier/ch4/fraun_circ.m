% ===== 书 4.5 节：fraun_circ —— 圆孔夫琅禾费衍射（Airy 图样） =====
% 书 4.5 节的主脚本（原书第 1-27 行），本文件忠实复现书图 4.3 并附数值检查。
% 参数（书 4.5 节）：w = 1 mm（圆孔半径，书注释写 "x half-width"）、
% lambda = 0.633 um（He-Ne 激光）、z = 50 m、L = 5 x 1.22*lambda*z/w ≈ 0.2 m、
% M = 250。辐照度用式 (4.35)：I2 = (w^2/(lambda z))^2 * jinc(w*rho/(lambda z))^2。
% 依赖：同文件夹的 jinc.m（书 4.5 节建的函数，J1(2*pi*x)/x）。
% MATLAB/Python 差异：
%   - e-3 与 10^-3 等价（书特意提醒：e 记号里不要用 ^ 符号！）；
%   - nthroot(I2,3) 是 3 次方根（拉出圆环），Python 写 I2**(1/3)；
%   - imagesc 默认 y 向下增，axis xy 翻回笛卡尔方向（Python 用 origin='lower'）。
clear; clc; close all;

% ---------- 书原代码（第 1-27 行） ----------
L=0.2;          %side length (m)
M=250;          %# samples
dx=L/M;         %sample interval
x=-L/2:dx:L/2-dx; y=x;  %coords
[X,Y]=meshgrid(x,y);

w=1e-3;         %x half-width
lambda=0.633e-6;%wavelength
z=50;           %prop distance
k=2*pi/lambda;  %wavenumber
lz=lambda*z;

%irradiance
I2=(w^2/lz)^2.*(jinc(w/lz*sqrt(X.^2+Y.^2))).^2;

figure(1)   %irradiance image
imagesc(x,y,nthroot(I2,3));
xlabel('x (m)'); ylabel('y (m)');
colormap('gray');
axis square;
axis xy;

figure(2)   %x-axis profile
plot(x,I2(M/2+1,:));
xlabel('x(m)'); ylabel('Irradiance');

% ---------- 数值检查（对照式 4.36/4.37 与菲涅耳数判据） ----------
fprintf('峰值 I2(0,0) = %.6e；解析 (w^2/(lam z))^2*pi^2 = %.6e\n', ...
        I2(M/2+1,M/2+1), (w^2/lz)^2*pi^2);
fprintf('第一零点解析 0.61*lam*z/w = %.5f m（式 4.36：2*pi*(w/lz)x = 1.22*pi）\n', ...
        0.61*lz/w);
fprintf('D_lobe = 1.22*lam*z/w = %.5f m；L = 5*D_lobe = %.4f m（书据此取 0.2 m）\n', ...
        1.22*lz/w, 5*1.22*lz/w);
fprintf('NF = w^2/(lam z) = %.4f << 1（判据 NF<0.1 -> z > 10*w^2/lam = %.2f m）\n', ...
        w^2/lz, 10*w^2/lam);
