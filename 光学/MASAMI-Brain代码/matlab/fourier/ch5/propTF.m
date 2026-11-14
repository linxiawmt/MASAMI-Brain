function[u2]=propTF(u1,L,lambda,z);
% ===== 书 5.1 节:菲涅耳传递函数(TF)传播器(逐行忠实 + 中文注释) =====
% 书 5.1 节 New M-file(命名为 "propTF")建立的函数(原书第 1-23 行):
%   1  function[u2]=propTF(u1,L,lambda,z);
%   11 [M,N]=size(u1);           %get input field array size
%   12 dx=L/M;                   %sample interval
%   13 k=2*pi/lambda;            %wavenumber
%   15 fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
%   16 [FX,FY]=meshgrid(fx,fx);
%   18 H=exp(-j*pi*lambda*z*(FX.^2+FY.^2));  %trans func
%   19 H=fftshift(H);            %shift trans func
%   20 U1=fft2(fftshift(u1));    %shift, fft src field
%   21 U2=H.*U1;                 %multiply
%   22 u2=ifftshift(ifft2(U2));  %inv fft, center obs field
%   23 end
% 要点(书 5.1 节注 a-g):
%   - (a) 第 11 行:size 只用 M,少传一个参数;
%   - (b) 第 16 行:fy 与 fx 相同,meshgrid 里 fx 用两次省一行;
%   - (c) 第 18 行:式 (5.2) 的 H 丢掉了 exp(jkz) 整体相位——它不影响
%     观察面的横向结构;
%   - (d)(e) 第 19/20 行:H 与 u1 都按"居中排列"存放,喂给 FFT 前先
%     fftshift 移位;
%   - (f)(g) 第 21/22 行:逐点相乘完成频域滤波,ifftshift 把 u2 移回居中
%     供显示。
% 记号提醒(书原话):小写 u 表示空域场、大写 U 表示傅里叶域量——与解析
%   式里大写表示空域场(如式 (5.1) 的 U1)不一致,"But what can you do?"
% MATLAB/Python 差异:Python 版在 python/fourier/ch5/ch5_prop.py,逐行对应;
%   numpy 索引 0-based,"零坐标样本"在 (M/2+1,M/2+1) 对应 (M//2,M//2)。
% ------------- 以下为书原代码(第 1-23 行) -------------
% propagation - transfer function approach
% assumes same x and y side lengths and
% uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field

[M,N]=size(u1);           %get input field array size
dx=L/M;                   %sample interval
k=2*pi/lambda;            %wavenumber

fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
[FX,FY]=meshgrid(fx,fx);

H=exp(-j*pi*lambda*z*(FX.^2+FY.^2));  %trans func
H=fftshift(H);            %shift trans func
U1=fft2(fftshift(u1));    %shift, fft src field
U2=H.*U1;                 %multiply
u2=ifftshift(ifft2(U2));  %inv fft, center obs field
end
