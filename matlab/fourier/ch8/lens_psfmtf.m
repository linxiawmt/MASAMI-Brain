% lens_psfmtt
% ===== 书 8.5 节:lens_psfmtf f/5 平凸透镜(Newport KPX094)的 PSF 与 MTF =====
% 书 8.5 节列出的脚本(原书第 1-68 行)。f/5 平凸透镜,ZEMAX 给出的 Seidel
% 系数(书表 8.2,单位:波长数 x 0.55 um)——注意书评注 (b) 写
% "coefficients taken from Table 8.1",应为表 8.2(勘误,见笔记)。
% 流程:式 (8.3)/(8.4) 相干传递函数 H = circ(f_r/f0)*exp(-jkW)
%       -> 式 (8.7) PSF = |IFFT{H}|^2 -> 式 (8.9) MTF = |FFT{PSF}| 归一化 DC
%       -> 与一维衍射限解析 MTF(式 7.31)对比。
% 采样判据(式 8.6):2f0 = 1/(lambda*fnum) = 3.6e5 1/m <= fN = 1/(2*du) = 5.1e5 1/m;
% 瞳直径占阵列 0.355(< 0.5 满足,书评注 c)。
% 轴上像点(u0=0,v0=0)改 (0,1) 即最大像高点:PSF 沿 v 拉长、MTF 分两个方向。
% MATLAB/Python 差异:Python 版 python/fourier/ch8/ch8_lens_psfmtf.py;
%   剖面用 h2(M/2+1,:)/( :,M/2+1)(1-based),numpy 为 h2[M//2,:]/h2[:,M//2];
%   nthroot(h2,2) 即 sqrt,imshow 对应 imagesc(...)+axis xy。
% ------------- 以下为书原代码(第 1-68 行) -------------
% lens_psfmtf
% f/5 plano-convex lens (Newport lens KPX094)
% 10 mm image height
% psf and mtf with Seidel aberrations
% aberration coefficients from ZEMAX
M=1024;             %sample #
L=1e-3;             %image plane side length
du=L/M;             %sample interval
u=-L/2:du:L/2-du; v=u; %coordinates

lambda=0.55e-6;     %wavelength
k=2*pi/lambda;      %wavenumber
Dxp=20e-3; wxp=Dxp/2; %exit pupil size
zxp=100e-3;         %exit pupil distance
fnum=zxp/(2*wxp);   %exit pupil f-number
lz=lambda*zxp;
twof0=1/(lambda*fnum);%incoh cutoff freq

u0=0; v0=0;         %normalized image coordinate

% aberration coefficients
wd=0*lambda;
w040=4.963*lambda;
w131=2.637*lambda;
w222=9.025*lambda;
w220=7.536*lambda;
w311=0.157*lambda;

fu=-1/(2*du):1/L:1/(2*du)-(1/L); %image freq coords
[Fu,Fv]=meshgrid(fu,fu);

% wavefront
W=seidel_5(u0,v0,-lz*Fu/wxp,-lz*Fv/wxp,...
    wd,w040,w131,w222,w220,w311);

% coherent transfer function
H=circ(sqrt(Fu.^2+Fv.^2)*lz/wxp).*exp(-j*k*W);
figure(1);
imagesc(u,v,angle(H)); axis xy; axis square
xlabel('u (m)'); ylabel('v (m)'); colormap('gray')

% point spread function
h2=abs(ifftshift(ifft2(fftshift(H)))).^2;

figure(2)       % psf image and profiles
imagesc(u,v,nthroot(h2,2)); axis xy; axis square
xlabel('u (m)'); ylabel('v (m)'); colormap('gray')
figure(3);
plot(u,h2(M/2+1,:)); xlabel('u (m)'); ylabel('PSF');
figure(4);
plot(u,h2(:,M/2+1)); xlabel('v (m)'); ylabel('PSF');

% MTF
MTF=fft2(fftshift(h2));
MTF=abs(MTF/MTF(1,1));  %normalize DC to 1
MTF=ifftshift(MTF);

% analytic MTF
MTF_an=(2/pi)*(acos(fu/twof0)-(fu/twof0)...
    .*sqrt(1-(fu/twof0).^2));
MTF_an=MTF_an.*rect(fu/(2*twof0)); %zero after cutoff

figure(5)       % MTF profiles
plot(fu,MTF(M/2+1,:),fu,MTF(:,M/2+1),':',...
    fu,MTF_an,'--');
axis([0 150000 0 1]);
legend('u MTF','v MTF','diff limit');
xlabel('f (cyc/m)'); ylabel('Modulation');
