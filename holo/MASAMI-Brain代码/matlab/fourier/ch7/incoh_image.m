% ===== 书 7.3.3 节:incoh_image 非相干成像例子(逐行忠实 + 中文注释) =====
% 书 7.3.3 节列出的脚本(原书第 1-48 行),实现式 (7.33):
%   Ii = F^-1{ OTF * F{Ig} },OTF = H 的零频归一化自相关(式 7.25/7.26)。
% 关键三行(第 22-23 行):OTF=ifft2(abs(fft2(fftshift(H))).^2);OTF=OTF/OTF(1,1);
%   ——FFT 自相关定理:|FT{H}|^2 的逆变换 = H 的(循环)自相关;
%   除以 OTF(1,1)(零频/零移位值)归一到 DC = 1(式 7.25);
%   OTF 故意留在 FFT 排列(不 fftshift),与第 32 行 fft2(fftshift(Ig))
%   的排列一致,直接逐点相乘(书评注 c)。
% 第 36 行:Ii=real(Ii); mask=Ii>=0; Ii=mask.*Ii; —— 数值噪声可能产生
%   微小的虚部/负值,imagesc/plot 不认,清零即可(书评注 d)。
% 图:figure(2) OTF 曲面(图 7.11)、figure(3) 模拟非相干像(图 7.12)、
%   figure(4) v=+0.2e-4 m 行剖面(图 7.13:组-1/单元3 刚可分辨,无振铃)。
% MATLAB/Python 差异:Python 版 ch7_incoh_image.py(ch7_imaging.py 的
%   otf_from_H/incoh_image 与第 22-23/32-36 行逐行对应);OTF(1,1) 是
%   1-based 的零频样本,numpy 写 OTF[0,0]。
% ------------- 以下为书原代码(第 1-48 行,第 3 行文件名有改动) -------------
%  incoh_image Incoherent Imaging Example

A=imread('USAF1951_py250','png'); %read image file
[M,N]=size(A);         %get image sample size
A=flipud(A);           %reverse row order
Ig=single(A);          %integer to floating
Ig=Ig/max(max(Ig));    %normalize ideal image

L=0.3e-3;              %image plane side length (m)
du=L/M;                %sample interval (m)
u=-L/2:du:L/2-du; v=u;

lambda=0.5*10^-6;      %wavelength
wxp=6.25e-3;           %exit pupil radius
zxp=125e-3;            %exit pupil distance
f0=wxp/(lambda*zxp);   %coherent cutoff

fu=-1/(2*du):1/L:1/(2*du)-(1/L); %freq coords
fv=fu;
[Fu,Fv]=meshgrid(fu,fv);
H=circ(sqrt(Fu.^2+Fv.^2)/f0);
OTF=ifft2(abs(fft2(fftshift(H))).^2);
OTF=OTF/OTF(1,1);

figure(2)              %check OTF
surf(fu,fv,fftshift(abs(OTF)))
camlight left; lighting phong
colormap('gray')
shading interp
xlabel('fu (cyc/m)'); ylabel('fv (cyc/m)');

Gg=fft2(fftshift(Ig)); %convolution
Gi=Gg.*OTF;
Ii=ifftshift(ifft2(Gi));
%remove residual imag parts, values <  0
Ii=real(Ii); mask=Ii>=0; Ii=mask.*Ii;

figure(3)              %image result
imagesc(u,v,nthroot(Ii,2));
colormap('gray'); xlabel('u (m)'); ylabel('v (m)');
axis square;
axis xy;

figure(4)              %horizontal image slice
vvalue=0.2e-4;         %select row (y value)
vindex=round(vvalue/du+(M/2+1)); %convert row index
plot(u,Ii(vindex,:),u,Ig(vindex,:),':');
xlabel('u (m)'); ylabel('Irradiance');
