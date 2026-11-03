% ===== 书 7.2.4 节:coh_image_rough 粗糙物体的相干成像(speckle)=====
% 书 7.2.4 节:把 coh_image 第 9 行换成随机相位版,重跑即得图 7.10:
%   ug=sqrt(Ig).*exp(j*2*pi*rand(M));   %ideal image(rough surface)
% rand(M) 给 [0,1) 均匀随机数,乘 2*pi 后每个样本获得 [0,2*pi) 独立相位
% ——"相对光波长极粗糙"的表面。相干叠加各点波前 -> speckle(相干斑)。
% 习题 7.5:
%   (a) 把 2*pi 换小(如 pi、0.4*pi):相位起伏幅度下降 -> 斑纹对比度
%       断崖下跌,测试板图案重新可辨(Python 实测均匀物体对比度
%       C=std/mean:2*pi -> 0.998(充分发展斑纹,理论 1)、pi -> 0.20、
%       0.4*pi -> 0.03);
%   (b) 改孔径(f/5 -> f/20):PSF 变宽 -> 斑纹颗粒变大,颗粒线度约
%       1.22*lam*(f/#)(Python 实测 3.6/6.0/9.6 um,网格 1.2 um)。
% MATLAB/Python 差异:rand 与 numpy default_rng 序列不同,复现图样
%   "统计一致、逐像素不同";Python 版 ch7_rough_speckle.py。
% ------------- 以下为书原代码 + 习题 7.5 的改动行 -------------
%  coh_image_rough Coherent Imaging Example (rough object)

A=imread('USAF1951_py250','png'); %read image file
[M,N]=size(A);         %get image sample size
A=flipud(A);           %reverse row order
Ig=single(A);          %integer to floating
Ig=Ig/max(max(Ig));    %normalize ideal image

ug=sqrt(Ig).*exp(j*2*pi*rand(M)); %ideal image(书 7.2.4 的唯一改动行)
% 习题 7.5(a):ug=sqrt(Ig).*exp(j*0.4*pi*rand(M));   %减弱粗糙度
L=0.3e-3;              %image plane side length (m)
du=L/M;                %sample interval (m)
u=-L/2:du:L/2-du; v=u;

lambda=0.5*10^-6;      %wavelength
wxp=6.25e-3;           %exit pupil radius(f/10;习题 7.5(b) 换 12.5e-3 -> f/5,
%                        3.125e-3 -> f/20)
zxp=125e-3;            %exit pupil distance
f0=wxp/(lambda*zxp);   %cutoff frequency

fu=-1/(2*du):1/L:1/(2*du)-(1/L); %freq coords
fv=fu;
[Fu,Fv]=meshgrid(fu,fv);
H=circ(sqrt(Fu.^2+Fv.^2)/f0);

H=fftshift(H);
Gg=fft2(fftshift(ug));
Gi=Gg.*H;
ui=ifftshift(ifft2(Gi));
Ii=(abs(ui)).^2;

figure(1)              %image result: speckle
imagesc(u,v,Ii);
colormap('gray'); xlabel('u (m)'); ylabel('v (m)');
axis square
axis xy
