function[u2]=propIR(u1,L,lambda,z);
% ===== 书 5.2 节:菲涅耳脉冲响应(IR)传播器(逐行忠实 + 中文注释) =====
% 书 5.2 节从 propTF 复制修改(New M-file 命名为 "propIR",原书第 1-23 行):
%   15 x=-L/2:dx:L/2-dx;         %spatial coords
%   16 [X,Y]=meshgrid(x,x);
%   18 h=1/(j*lambda*z)*exp(j*k/(2*z)*(X.^2+Y.^2)); %impulse
%   19 H=fft2(fftshift(h))*dx^2; %create trans func
% 其余与 propTF 相同。
% 要点(书 5.2 节注 a/b):
%   - (a) 第 18 行:实现式 (5.4) 的脉冲响应 h,再用 FFT 得到 H;
%   - (b) 第 19 行:乘 dx^2 是关键——U1 的 FFT 与 U2 的 IFFT 定标互相
%     抵消,唯有 FFT{h} 这一步需要自己的定标因子 dx^2 才能回到连续量纲。
% 解析上与 propTF 等价(式 (5.1) 与 (5.3) 是同一运算),但离散采样后
% 两个传播器的行为可以不同(书 5.4 节的主题):
%   - Δx > λz/L 时 h 的 chirp 相位欠采样 -> 观察面出现周期副本;
%   - Δx < λz/L 时 h 过采样,FFT 兼做"窗口滤波",压掉 L/(2λz) 以外的
%     源频谱 -> 长距离下更平滑但不总准确。
% MATLAB/Python 差异:Python 版在 python/fourier/ch5/ch5_prop.py;
%   MATLAB 偶数阵列 fftshift 与 ifftshift 相同(书全用偶数 M)。
% ------------- 以下为书原代码(第 1-23 行) -------------
% propagation - impulse response approach
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

x=-L/2:dx:L/2-dx;         %spatial coords
[X,Y]=meshgrid(x,x);

h=1/(j*lambda*z)*exp(j*k/(2*z)*(X.^2+Y.^2)); %impulse
H=fft2(fftshift(h))*dx^2; %create trans func
U1=fft2(fftshift(u1));    %shift, fft src field
U2=H.*U1;                 %multiply
u2=ifftshift(ifft2(U2));  %inv fft, center obs field
end
