function[u2,L2]=propFF(u1,L1,lambda,z);
% ===== 书 5.5 节:夫琅禾费传播器(逐行忠实 + 中文注释) =====
% 书 5.5 节 New M-file(命名为 "propFF")建立的函数(原书第 1-22 行):
%   11 [M,N]=size(u1);           %get input field array size
%   12 dx1=L1/M;                 %source sample interval
%   13 k=2*pi/lambda;            %wavenumber
%   15 L2=lambda*z/dx1;          %obs sidelength
%   16 dx2=lambda*z/L1;          %obs sample interval
%   17 x2=-L2/2:dx2:L2/2-dx2;    %obs coords
%   18 [X2,Y2]=meshgrid(x2,x2);
%   20 c=1/(j*lambda*z)*exp(j*k/(2*z)*(X2.^2+Y2.^2));
%   21 u2=c.*ifftshift(fft2(fftshift(u1)))*dx1^2;
%   22 end
% 要点:
%   - 观察面参数由式 (5.18) 决定:L2 = lambda*z/dx1、dx2 = lambda*z/L1
%     ——源面与观察面边长一般不同!函数额外输出 L2,省得主脚本重算
%     (书原话)。
%   - 若临界采样(dx1 = lambda*z/L1),则 L2 = L1;否则两边长不同。
%   - 第 21 行:fftshift(u1) 移位 -> fft2 -> ifftshift 移回居中,
%     再乘前导因子 c 与定标 dx1^2(连续 FT 的 DFT 近似,书 2.3 节)。
%   - 只关心辐照度时前导相位因子 e^{jkz} 与 chirp c 都可以丢(书 5.5 节
%     末尾);但要注意 chirp 的采样条件(式 (5.23))。
% MATLAB/Python 差异:Python 版在 python/fourier/ch5/ch5_prop.py 的
%   propFF(返回 (u2, L2) 元组)。
% ------------- 以下为书原代码(第 1-22 行) -------------
% propagation - Fraunhofer pattern
% assumes uniform sampling
% u1 - source plane field
% L1 - source plane side length
% lambda - wavelength
% z - propagation distance
% L2 - observation plane side length
% u2 - observation plane field

[M,N]=size(u1);           %get input field array size
dx1=L1/M;                 %source sample interval
k=2*pi/lambda;            %wavenumber

L2=lambda*z/dx1;          %obs sidelength
dx2=lambda*z/L1;          %obs sample interval
x2=-L2/2:dx2:L2/2-dx2;    %obs coords
[X2,Y2]=meshgrid(x2,x2);

c=1/(j*lambda*z)*exp(j*k/(2*z)*(X2.^2+Y2.^2));
u2=c.*ifftshift(fft2(fftshift(u1)))*dx1^2;
end
