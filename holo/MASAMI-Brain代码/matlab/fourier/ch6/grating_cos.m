% ===== 书 6.4.1 节:grating_cos 余弦振幅光栅例子(逐行忠实 + 中文注释) =====
% 参数:lambda=0.5 um、f=0.5 m、P=1e-4 m(光栅周期)、D1=1.02e-3(光栅边长,
%   跨 rect 函数 D1/dx1=51 个样本——奇数,与 rect 定义一致)、
%   L1=1e-2 m、M=500。跨每个余弦周期的样本数 M*P/L1=5(式 6.21,"okay")。
% 图 6.10:(a) 源面辐照度(显示分辨率看不出周期结构);(b) 夫琅禾费图样
%   (三次方根拉伸);(c) x 剖面——中央零阶 + ±1 一阶峰(位置 ±lambda*f/P=±2.5 mm)。
% 解析对照:式 (6.22)/(6.23) 的 FT 与式 (6.24) 的夫琅禾费场(丢掉前导复指数,
%   只算辐照度,"but don't forget 1/lambda*z",这里 z 就是 f)。
% 波长分离试验(图 6.11):换 lambda=0.6e-6 重跑,一阶峰移到 ±3.0 mm;
%   两条曲线的 x2 网格不同,不能用同一组 x2 画图(书特别提醒)!
% MATLAB/Python 差异:Python 版 ch6_gratings.py 的 grating_cos 部分;
%   sinc 为信号处理工具箱函数,无工具箱可自定义 sinc_v=@(x) sin(pi*x)./(pi*x)。
% ------------- 以下为书原代码(第 1-21 行 + 书给的分析代码) -------------
% grating_cos diffraction grating example

lambda=0.5e-6;  %wavelength
f=0.5;          %focal distance
P=1e-4;         %grating period
D1=1.02e-3;     %grating side length

L1=1e-2;        %array side length
M=500;          %# samples
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; %source coords
[X1,Y1]=meshgrid(x1,x1);

% Grating field and irradiance
u1=1/2*(1-cos(2*pi*X1/P)).*rect(X1/D1).*rect(Y1/D1);

% Fraunhofer pattern
[u2,L2]=propFF(u1,L1,lambda,f);
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; y2=x2;%obs coords
I2=abs(u2).^2;
%
figure(1)             %图 6.10(a):源面辐照度
imagesc(x1,y1,abs(u1).^2); axis square; axis xy; colormap('gray');
xlabel('x (m)'); ylabel('y (m)');
figure(2)             %图 6.10(b):夫琅禾费图样(三次方根拉伸)
imagesc(x2,y2,nthroot(I2,3)); axis square; axis xy; colormap('gray');
xlabel('x (m)'); ylabel('y (m)');
figure(3)             %图 6.10(c):x 剖面
plot(x2,I2(M/2+1,:)); xlabel('x (m)'); ylabel('Irradiance');

% analytic(书 6.4.1 节第二段代码,式 (6.24) 的实现)
[X2,Y2]=meshgrid(x2,y2);
lf=lambda*f;
u2a=(1/lf)*D1^2/2*sinc(D1/lf*Y2)...
    .*(sinc(D1/lf*X2)-1/2*sinc(D1/lf*(X2+lf/P))...
    -1/2*sinc(D1/lf*(X2-lf/P)));
I2a=abs(u2a).^2;
figure(4)             %数值(实线)vs 解析(虚线)
plot(x2,I2(M/2+1,:)); hold on;
plot(x2,I2a(M/2+1,:),'--'); hold off;
xlabel('x (m)'); ylabel('Irradiance');
legend('discrete','analytic');
