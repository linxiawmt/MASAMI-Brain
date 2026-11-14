% ===== 书 7.2.3 节:coh_image 相干成像例子(逐行忠实 + 中文注释) =====
% 书 7.2.3 节列出的脚本(原书第 1-51 行),实现式 (7.16):
%   Ui = F^-1{ H(fU,fV) * F{Ug} },H = circ(sqrt(Fu^2+Fv^2)/f0)(式 7.13)。
% 参数:f/10 薄透镜(wxp=6.25 mm、zxp=125 mm)、lambda=0.5 um、
%   L=0.3 mm、M=250 -> du=1.2 um <= lam*(f/#)/2 = 2.5 um(式 7.18)。
% 频率盘点:fN = 4.17e5、2f0 = 2e5、f0 = 1e5 cyc/m,判据 2f0 <= fN(式 7.17)。
% 原书读 USAF1951B250.png(书网站下载,已失效);本目录先跑
% make_usaf_chart.m 生成同规格 USAF1951_py250.png,这里改读它。
% 图:figure(1) 理想图像(图 7.5)、figure(2) 相干传递函数(图 7.6)、
%   figure(3) 模拟相干像 sqrt(Ii) 拉伸(图 7.7)、figure(4) 行剖面(图 7.8,
%   v=-0.8e-4 m:组-2/单元6 条纹周期约 1e-5 m,恰为相干截止 f0)。
% MATLAB/Python 差异:Python 版 ch7_coh_image.py;
%   - 第 30 行 H.*.99 是 MATLAB 7.1 画单位高圆面的显示权宜(surf 画不动
%     单位高度),新版 MATLAB 可去掉;
%   - 第 34 行 xlabel('fv')/ylabel('fu') 与 surf(fu,fv,...) 的轴映射相反
%     (surf 第一参数 fu 应配 x 轴)——圆瞳下无影响,矩形瞳(习题 7.3)
%     必须改正,本文件已按正确映射写;
%   - 第 49 行 vindex=round(vvalue/du+(M/2+1)):1-based;numpy 是
%     M//2 + round(vvalue/du)。
% ------------- 以下为书原代码(第 1-51 行,第 3 行文件名有改动) -------------
%  coh_image Coherent Imaging Example

A=imread('USAF1951_py250','png'); %read image file
[M,N]=size(A);         %get image sample size
A=flipud(A);           %reverse row order
Ig=single(A);          %integer to floating
Ig=Ig/max(max(Ig));    %normalize ideal image

ug=sqrt(Ig);           %ideal image field
L=0.3e-3;              %image plane side length (m)
du=L/M;                %sample interval (m)
u=-L/2:du:L/2-du; v=u;

figure(1)              %check ideal image
imagesc(u,v,Ig);
colormap('gray');  xlabel('u (m)'); ylabel('v (m)');
axis square
axis xy

lambda=0.5*10^-6;      %wavelength
wxp=6.25e-3;           %exit pupil radius
zxp=125e-3;            %exit pupil distance
f0=wxp/(lambda*zxp);   %cutoff frequency

fu=-1/(2*du):1/L:1/(2*du)-(1/L); %freq coords
fv=fu;
[Fu,Fv]=meshgrid(fu,fv);
H=circ(sqrt(Fu.^2+Fv.^2)/f0);

figure(2)              %check H
surf(fu,fv,H.*.99)
camlight left; lighting phong
colormap('gray')
shading interp
xlabel('fu (cyc/m)'); ylabel('fv (cyc/m)');   %书原文两条标签互换,圆瞳无影响

H=fftshift(H);
Gg=fft2(fftshift(ug));
Gi=Gg.*H;
ui=ifftshift(ifft2(Gi));
Ii=(abs(ui)).^2;

figure(3)              %image result
imagesc(u,v,nthroot(Ii,2));
colormap('gray'); xlabel('u (m)'); ylabel('v (m)');
axis square
axis xy

figure(4)              %horizontal image slice
vvalue=-0.8e-4;       %select row (y value)
vindex=round(vvalue/du+(M/2+1)); %convert row index
plot(u,Ii(vindex,:),u,Ig(vindex,:),':');
xlabel('u (m)'); ylabel('Irradiance');
