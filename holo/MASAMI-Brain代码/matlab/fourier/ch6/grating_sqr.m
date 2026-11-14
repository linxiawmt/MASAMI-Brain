% ===== 书 6.4.2 节:grating_sqr 方波振幅光栅例子(逐行忠实 + 中文注释) =====
% 透过率(式 6.25):t_A = [rect(x1/(P/2)) * (1/P)comb(x1/P)] rect(x1/D1) rect(y1/D1)
%   ——方波 = 宽 P/2 的 rect 与周期 P 的梳做一维卷积,再乘 D1 x D1 窗。
% 参数:lambda=0.5 um、f=0.5 m、P=1e-4 m、D1=1e-3 m、L1=1e-2 m、M=1000。
%   dx1=1e-5 m,P/dx1=10(整数!ucomb 落位的前提),P/2=5(奇数个样本跨
%   每个 rect,与 rect 定义一致)。
% 卷积定标(书原话):结果不需要再乘定标因子——(a) ifft 已为一次 fft 正确
%   定标(例如 fr);(b) 每个单位样本 delta 表现得像面积 1/dx1,恰好抵消
%   第二个 fft 带来的 dx1。repmat 把一维 ux 快速铺满二维 u1 的每一行。
% 图 6.13:与余弦光栅类似但多出高阶峰(奇数阶 n=±1,±3,±5;偶数阶被
%   sinc(n/2) 的零点压掉)。解析对照:式 (6.27),n=-5:5 共 11 项足够。
% MATLAB/Python 差异:Python 版 ch6_gratings.py 的 grating_sqr 部分;
%   ucomb 来自附录 C.5(本目录提供 ucomb.m);书末尾"Implement this code
%   at the end of sqr_grating"中的 sqr_grating 是 grating_sqr 之笔误。
% ------------- 以下为书原代码(第 1-19 行 + 书给的解析代码) -------------
% grating_sqr diffraction grating example

lambda=0.5e-6;  %wavelength
f=0.5;          %focal distance
P=1e-4;         %grating period
D1=1e-3;        %grating side length

L1=1e-2;        %array side length
M=1000;         %# samples
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; %source coords
[X1,Y1]=meshgrid(x1,x1);

% construct grating field
fc=fft(fftshift(ucomb(x1/P)));
fr=fft(fftshift(rect(x1/(P/2))));
ux=ifftshift(ifft(fc.*fr)); %1D conv rect & comb
u1=repmat(ux,M,1);          %replicate to 2D
u1=u1.*rect(X1/D1).*rect(Y1/D1); %set size

... (then apply the Fraunhofer propagator)   %书原文此处省略——实际脚本:
[u2,L2]=propFF(u1,L1,lambda,f);
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; y2=x2;%obs coords
I2=abs(u2).^2;
%
figure(1)                    %图 6.13:x 剖面
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');

% analytic(书 6.4.2 节代码,式 (6.27) 的 y=0 切片,n=-5:5)
lf=lambda*f;
u2a=zeros(1,M);
for n=-5:5
    ut=sinc(n/2)*sinc(D1/lf*(x2-n*lf/P));
    u2a=u2a+ut;
end
u2a=D1^2/(2*lf)*u2a;
I2a=abs(u2a).^2;
figure(2)                    %数值(实线)vs 解析(虚线)
plot(x2,I2(M/2+1,:)); hold on;
plot(x2,I2a,'--'); hold off;
xlabel('x (m)'); ylabel('Irradiance');
legend('discrete','analytic');
