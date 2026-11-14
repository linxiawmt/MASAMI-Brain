% ===== 书 6.4.3 节:grating_sqr1D 一维光栅模型(逐行忠实 + 中文注释) =====
% 光栅没有 y 方向周期性时,可以只算 y=0 的一维切片:y 相关项删掉、
% meshgrid 不再需要、repmat 不再需要、二维 propFF 换成一维夫琅禾费。
% 好处:同样的内存能放下更长的向量——更大的总宽度、跨光栅更多周期。
% 参数:物理参数与 grating_sqr 完全相同,M=2200(书指定值:
%   P/dx1 = M*P/L1 = 22,整数;每个 rect 跨 P/2 = 11 个样本,奇数)。
% 一维夫琅禾费(书给的代码):前导因子 1/sqrt(lambda*f)——二维是
%   1/(lambda*f),一维开平方根;丢掉 e^{jkf} 与 chirp(只关心辐照度)。
%   L2 = lambda*f/dx1 = 0.055 m,比二维(M=1000)的 0.05 m 更宽。
% 图 6.14:(a) 全向量视图;(b) 局部放大——阶的位置与相对大小与二维
%   结果一致;幅度量级不同(一维模型不计第二维的功率)。
% 注意:书正文说"The expanded view in Fig. 6.11(b)"是笔误,应为图 6.14(b)。
% MATLAB/Python 差异:Python 版 ch6_gratings.py 的 grating_sqr1D 部分。
% ------------- 以下为书脚本要点(一维改写) -------------
% grating_sqr1D 1D square-wave grating model
lambda=0.5e-6;  %wavelength
f=0.5;          %focal distance
P=1e-4;         %grating period
D1=1e-3;        %grating side length

L1=1e-2;        %array side length
M=2200;         %# samples
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; %source coords(一维向量)

% construct grating field(一维:无 meshgrid / repmat)
fc=fft(fftshift(ucomb(x1/P)));
fr=fft(fftshift(rect(x1/(P/2))));
ux=ifftshift(ifft(fc.*fr)); %1D conv rect & comb
u1=ux.*rect(x1/D1);         %set size(不再乘 rect(Y1/D1))

% 1D Fraunhofer pattern(书给的三行)
lf=lambda*f;
u2=sqrt(1/lf)*ifftshift(fft(fftshift(u1)))*dx1;
L2=lf/dx1;
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2;      %obs coords
I2=abs(u2).^2;
%
figure(1)                   %图 6.14(a):全向量视图
plot(x2,I2);
xlabel('x (m)'); ylabel('Irradiance');
%
figure(2)                   %图 6.14(b):局部放大
plot(x2,I2); xlim([-0.01 0.01]);
xlabel('x (m)'); ylabel('Irradiance');
