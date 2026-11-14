function [Dphi rgrid] = str_fcn2_ft(phz, delta)
% 二维相位结构函数（书 Listing 3.7 的复现；ch3 分册 PDF 缺失，
% 按 9.5.5 节的用法重写，与 Python 版 python/waveprop/ch9 数值一致）
%
% 按定义计算无偏结构函数：
%   D(rho) = <(phi(r + rho) - phi(r))^2> = [S1 + S2 - 2*C]/Npairs
%   S1(rho) = sum_{r in overlap} phi(r)^2，S2(rho) = sum phi(r+rho)^2，
%   C(rho)  = sum phi(r)*phi(r+rho)，Npairs = 重叠像素数（方形支撑）。
% 实现：相位与窗口填零到 2N x 2N 后做循环相关——支撑在 [0,N)^2、
% |rho| < N 时循环相关 == 线性相关（无卷绕污染），覆盖原书图 9.3
% 横轴 12*r0（约 0.75*N*delta）的大间隔需求。
%
% 注意：不要用 D = 2*(var(phi) - cov(rho)) 的"方差-协方差"捷径——
% 次谐波屏含波长超过屏幕的大幅低频模，全局方差不等于逐对平均的
% 均方，小间隔处会高估数倍。
%
% 用法：[Dphi rgrid] = str_fcn2_ft(phz, delta)
%   Dphi:  N×N 结构函数，零间隔在 (1,1)，只含 rho_x, rho_y >= 0 象限
%   rgrid: 与 Dphi 同尺寸的径向间隔网格 [m]

N  = size(phz, 1);
N2 = 2*N;

% 填零到 2N（FFT 序，支撑放左上角 [1:N, 1:N]）
P0  = zeros(N2); P0(1:N, 1:N)   = phz;
W0  = zeros(N2); W0(1:N, 1:N)   = ones(N);
P20 = zeros(N2); P20(1:N, 1:N)  = phz.^2;

FP  = fft2(P0);
FW  = fft2(W0);
FP2 = fft2(P20);

% wrap 序相关：Cfull(m) = sum_j A(j) B(j+m)
C   = real(ifft2(conj(FP)  .* FP));    % sum phi(j) phi(j+m)
S1  = real(ifft2(conj(FP2) .* FW));    % sum phi(j)^2 W(j+m)
S2  = real(ifft2(FP2       .* conj(FW)));  % sum W(j) phi(j+m)^2
Np  = real(ifft2(conj(FW)  .* FW));    % 重叠像素数

Dfull = (S1 + S2 - 2*C) ./ max(Np, 1);

% fftshift 后零间隔落在 (N+1, N+1)（1-based）；取 rho >= 0 象限
Dc   = fftshift(Dfull);
Dphi = Dc(N+1 : 2*N, N+1 : 2*N);
[ii jj] = ndgrid(0:N-1, 0:N-1);
rgrid = delta * sqrt(ii.^2 + jj.^2);
end
