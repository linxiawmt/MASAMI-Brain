% ===== 书 1.1 例 1.2 · Table 1.4 fft2Dbitmap_image.m =====
%fft2Dbitmap_image.m
%Simulation of Fourier transformation of bitmap images
% 对一幅 256x256 位图做二维傅里叶变换并显示幅度谱。
% 需要事先用画图板做一幅 256x256 的 triangle.bmp（黑色背景、白色三角形）放在同目录。
clear

I=imread('triangle.bmp','bmp'); %Input bitmap image   % 读入位图
I=I(:,:,1);                                            % 取一个通道（灰度图）
figure(1) %displaying input                            % 显示输入图像
colormap(gray(255));
image(I)
axis off
FI=fft2(I);                                            % 二维 FFT
FI=fftshift(FI);                                       % 零频搬到中心
max1=max(FI);
max2=max(max1);
scale=1.0/max2;                                        % 归一化到峰值 1
FI=FI.*scale;
figure(2) %Gray scale image of the absolute value of transform
colormap(gray(255));
image(10*(abs(256*FI)));                               % 幅度放大后按灰度显示（星芒状频谱）
axis off
