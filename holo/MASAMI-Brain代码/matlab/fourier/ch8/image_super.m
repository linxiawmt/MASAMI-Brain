% image_super superposition image
% ===== 书 8.7.2 节:image_super 空间变叠加成像仿真 =====
% 书 8.7.2 节列出的脚本(原书第 1-60 行)。对理想像 Ig 的每个样本 (u0,v0):
% 生成该像点的 PSF(式 8.7)-> 体积归一到 1(评注 e:等价式 7.25 的 OTF
% 归一化)-> circshift 平移到像点位置(评注 f)-> 乘 Ig 后累加(式 8.17 的
% 离散实现,评注 g)。
% 注意:原书读入随书代码分发的 USAF1951B250.png;Python 侧笔记用合成
% USAF 风格靶替代(ch8_aberrations.make_usaf_target)。
% 评注 (c)/(d):外层循环索引 n 是行(= y/v 方向),内层 m 是列(u 方向)
% ——本书少有的"用索引推物理参数"的场合,因为第 53 行要逐元素寻址 Ig。
% 评注 (h):M=250 共 62,500 个 PSF,书作者笔记本约 1 小时!对比:卷积
% 只要一个 PSF + 几次 FFT,瞬间完成——空间变系统没有这个奢侈。
% MATLAB/Python 差异:Python 版 python/fourier/ch8/ch8_image_super.py
%   (可传 M 参数快速验证;M=100 时式 8.6 不再满足,见笔记)。
% ------------- 以下为书原代码(第 1-60 行) -------------
% image_super superposition image

A=imread('USAF1951B250','png');

[M,N]=size(A); A=flipud(A);
Ig=single(A); Ig=Ig/max(max(Ig));

L=1e-3;                %image plane side length
du=L/M;                %sample interval
u=-L/2:du:L/2-du; v=u; %coordinates
fN=1/(2*du)            %Nyquist frequency

lambda=0.5*10^-6;      %wavelength
k=2*pi/lambda;         %wavenumber
wxp=2.5e-3;            %exit pupil radius
zxp=100e-3;            %exit pupil distance
fnum=zxp/(2*wxp)       %exit pupil f-number

twof0=1/(lambda*fnum)  %inc cutoff freq
fN=1/(2*du)            %Nyquist frequency

% aberration coefficients
wd=0*lambda;
w040=0.5*lambda;
w131=1*lambda;
w222=1.5*lambda;
w220=0*lambda;
w311=0*lambda;

% frequency coordinates
fu=-1/(2*du):1/L:1/(2*du)-(1/L);
fu=fftshift(fu); %shift cords, avoid shift H in loop
[Fu,Fv]=meshgrid(fu,fu);

I=zeros(M);
% loop through image plane positions
for n=1:M
   v0=(n-(M/2+1))/(M/2) %norm v image coord
   for m=1:M
       u0=(m-(M/2+1))/(M/2); %norm u image coord
       % wavefront
       W=seidel_5(u0,v0,-2*lambda*fnum*Fu,...
           -2*lambda*fnum*Fv,...
           wd,w040,w131,w222,w220,w311);
       % coherent transfer function
       H=circ(sqrt(Fu.^2+Fv.^2)*2*lambda*fnum)...
           .*exp(-j*k*W);
       % create PSF - normalize volume to 1
       h2=abs(ifftshift(ifft2(H))).^2;
       h2=h2/(sum(sum(h2)));
       % shift h2 to image plane position
       h2=circshift(h2,[n-(M/2+1),m-(M/2+1)]);
       % superposition
       I=Ig(n,m)*h2+I;
   end
end

figure(1)
imagesc(u,v,nthroot(I,3));
colormap('gray'); axis square; axis xy
xlabel('x (m)'); ylabel('y (m)');
