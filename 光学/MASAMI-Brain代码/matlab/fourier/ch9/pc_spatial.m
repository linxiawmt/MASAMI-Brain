% pc_spatial partial spatial coherence example
% ===== 书 9.2.3 节:pc_spatial 部分空间相干(随机相位屏 + 系综平均) =====
% 书 9.2.3 节列出的脚本(原书第 1-61 行),逐行忠实。与 pc_temp 同一光学布局
% 但去掉路径延迟:l_cr=8mm 的部分空间相干让条纹对比度下降(书图 9.11)。
% 参数链(书 21-22 行):sigma_f = 2.5*l_cr = 20mm,sigma_r = sqrt(4*pi*sigma_f^4
%   /l_cr^2) = 177.2mm(式 9.23 反解);式 9.21 的 sigma_r^2/(2*pi*sigma_f^2)
%   = 12.5 >> 1。
% 相位屏(书 39-40 行,式 9.25):fie = ifft2(F.*(randn+j*randn))*sigma_r/dfx1
%   *M^2*dfx1^2 —— sigma_r/dfx1 来自式 9.25 的方差归一(均匀采样下
%   (dfx*dfy)^(1/2) = dfx1;书评注 (e) 印刷版漏了 sigma_r),M^2 抵消 ifft2
%   内置的 1/M^2,dfx1^2 补上逆变换积分的面积元;
%   实部/虚部是两块独立屏("一次 FFT 两块屏")。
% 评注 (a)/(c):第 9/26 行对坐标数组本身 fftshift -> 场与滤波器 F 天生在
%   FFT 排列,循环内不用再移位。
% 解析对照(式 9.27):书印刷版花括号内是 w^2,应为 w^4(同式 9.10 的勘误;
%   按 w^4 中心值 1253,与书图 9.11(b) 峰值 ~1250 吻合)。
% 注意(数值验证发现,见笔记"数值验证"):L1=50mm 时 dfx1=20 1/m,而滤波器
%   F 的通带 1/(pi*sigma_f) = 15.9 1/m —— 通带不足一个频率网格,屏的结构
%   由 DC 与 ±20 cyc/m 少数样本决定,相关长度被拉长,条纹可见度实测 ~0.8
%   (解析 |mu_12| = exp(-(5/8)^2) = 0.68)。书 9.2.3 节描述的 "passband
%   defined by a few samples" 失效模式在 sigma_f < L/2 时也会发生;补救:
%   调小 21 行的 2.5 倍率(~1.2-1.5)或加大 L1 与 M 让通带被充分采样。
% 依赖:circ.m(书附录 C.1;ch7 目录有副本)。
% MATLAB/Python 差异:Python 版 python/fourier/ch9/ch9_spatial.py,逐行对应;
%   randn(M)+j*randn(M) 对应 default_rng().standard_normal((M,M))+1j*(...)
%   (Python 侧固定随机种子便于复现,书不设种子);剖面 I2(M/2+1,:)
%   对应 I2[M//2,:]。
% ------------- 以下为书原代码(第 1-61 行) -------------
% pc_spatial partial spatial coherence example

lambda=650e-9;      %center wavelength (m)

L1=50e-3;           %source plane side length
M=250;              %# samples (even)
dx1=L1/M;           %sample interval
x1=-L1/2:dx1:L1/2-dx1; %source coords
x1=fftshift(x1);    %shift x coord
[X1,Y1]=meshgrid(x1,x1);

% beam parameters
w=1e-3;             %radius
dels=5e-3;          %transverse separation
f=0.25;             %Fraunhofer focal distance
lf=lambda*f;

% partial spatial coherence screen parameters
N=100;              %number of screens (even)
Lcr=8e-3;           %spatial correlation length
sigma_f=2.5*Lcr;    %Gaussian filter parameter
sigma_r=sqrt(4*pi*sigma_f^4/Lcr^2); %random std

dfx1=1/L1;
fx1=-1/(2*dx1):dfx1:1/(2*dx1)-dfx1;
fx1=fftshift(fx1);
[FX1,FY1]=meshgrid(fx1,fx1);

% source field
u1=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)...
    +circ(sqrt((X1+dels/2).^2+Y1.^2)/w);
% filter spectrum
F=exp(-pi^2*sigma_f^2*(FX1.^2+FY1.^2));

% loop through screens
I2=zeros(M);
for n=1:N/2
   % make 2 random screens
   fie=(ifft2(F.*(randn(M)+j*randn(M)))...
       *sigma_r/dfx1)*M^2*dfx1^2;
   % Fraunhofer pattern applying screen 1
   u2=1/lf*(fft2(u1.*exp(j*real(fie))))*dx1^2;
   I2=I2+abs(u2).^2;
   % Fraunhofer pattern applying screen 2
   u2=1/lf*(fft2(u1.*exp(j*imag(fie))))*dx1^2;
   I2=I2+abs(u2).^2;
end

I2=ifftshift(I2)/N;    %normalize & center irradiance
x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf; %obs coords
y2=x2;

figure(1)              %irradiance image
imagesc(x2,y2,I2);
xlabel('x (m)'); ylabel('y (m)');
axis square; axis xy;
colormap('gray');

figure(2)              %irradiance slice
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
