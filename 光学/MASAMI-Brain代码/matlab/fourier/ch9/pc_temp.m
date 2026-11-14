% pc_temp partial temporal coherence example
% ===== 书 9.1.3 节:pc_temp 部分时间相干(双光束杨氏干涉 + 高斯谱线型) =====
% 书 9.1.3 节列出的脚本(原书第 1-57 行),逐行忠实。参数:lambda0=650nm、
% dnu(FWHM)=2GHz、w=1mm、dels=5mm、f=0.25m、deld=5cm -> tau = deld/c,
% |gamma(tau)| = exp[-(pi*dnu*tau/(2*sqrt(ln2)))^2] = 0.6733(书式 9.4/9.11)。
% 流程:对 N=51 条谱线(跨 4*dnu)各生成"带延迟相位的双圆孔"场(书 37-39 行,
%   延迟 exp(j*k*deld) 乘在 x=-dels/2 一束上)-> Fraunhofer 图样(书 41 行,
%   忽略前端相位,中心波长定标)-> 按 S(nu)*dnu 加权累加(书 43 行,式 9.9)。
% 评注 (a):51 条线跨 4*dnu。注意书评注说"This is the span displayed in
%   Fig. 9.2"——4*dnu 的采样谱线跨度画在图 9.3 里,图 9.2 是 |gamma(tau)|
%   曲线(tau 轴),图号笔误(见笔记勘误)。
% 评注 (c):第 19 行对"坐标数组本身"做 fftshift -> 场数组天生 DC 在索引 1,
%   循环内的 fft2 不用再移位(省时;与 ch8 psf_map 同一招)。
% 评注 (f):Fraunhofer 定标用中心波长 lambda0——谱宽对 1/lf 与 x2 的影响
%   可忽略,但光程差项 exp(j*k*deld) 里的 k 必须逐线更新(时间相干的全部
%   来源就在这个相位;若统一用 k0 会得到完全相干的结果)。
% 解析对照(式 9.10):书印刷版花括号内是 w^2,应为 w^4 —— 数值反证见
%   python/fourier/ch9/ch9_temporal.py(按 w^2 中心值 7.5e8,按 w^4 得
%   ~1250,与书图 9.5(b) 峰值 ~1200 吻合)。
% 依赖:circ.m(书附录 C.1;ch7 目录有副本)。
% MATLAB/Python 差异:Python 版 python/fourier/ch9/ch9_temporal.py,逐行对应;
%   circ(sqrt(...)/w) 对应 base.circle(X,Y,w)(Goodman 半径约定,见 ch1 勘误);
%   for n=1:N 对应 range(1,N+1);(N+1)/2 两种语言里都是 26.0;
%   剖面 I2(M/2+1,:) (1-based) 对应 I2[M//2,:] (0-based)。
% ------------- 以下为书原代码(第 1-57 行) -------------
% pc_temp partial temporal coherence example

lambda0=650e-9;     %center wavelength (m)
c=3e8;              %speed of light
k0=2*pi/lambda0;    %center wavenumber
nu0=c/lambda0;      %center frequency

% Gaussian lineshape parameters
N=51;               %number of components (odd)
delnu=2e9;          %spectral density FWHM (Hz)
b=delnu/(2*sqrt(log(2))); %FWHM scaling
dnu=4*delnu/N;      %freq interval

% source plane parameters
L1=50e-3;           %source plane side length
M=250;              %# samples (even)
dx1=L1/M;           %sample interval
x1=-L1/2:dx1:L1/2-dx1; %source coords
x1=fftshift(x1);     %shift x coord
[X1,Y1]=meshgrid(x1,x1);

% beam parameters
w=1e-3;             %radius
dels=5e-3;          %transverse separation
deld=5e-2;          %delay distance
f=0.25;             %focal dist for Fraunhofer
lf=lambda0*f;

% loop through lines
I2=zeros(M);
for n=1:N
   % spectral density function
   nu=(n-(N+1)/2)*dnu+nu0;
   S=1/(sqrt(pi)*b)*exp(-(nu-nu0)^2/b^2);
   k=2*pi*nu/c;
   % source
   u=circ(sqrt((X1-dels/2).^2+Y1.^2)/w)...
       +circ(sqrt((X1+dels/2).^2+Y1.^2)/w)...
       *exp(j*k*deld);
   % Fraunhofer pattern
   u2=1/lf*(fft2(u))*dx1^2;
   % weighted irradiance and sum
   I2=I2+S*(abs(u2).^2)*dnu;
end

I2=ifftshift(I2); %normalize/center irradiance
x2=(-1/(2*dx1):1/L1:1/(2*dx1)-1/L1)*lf; %obs coords
y2=x2;

figure(1)              %irradiance image
imagesc(x2,y2,I2);
xlabel('x (m)'); ylabel('y (m)');
axis square; axis xy; colormap('gray');

figure(2)              %irradiance profile
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
