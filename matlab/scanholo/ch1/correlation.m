% ===== 书 1.2 例 1.6 · Table 1.5 correlation.m =====
%correlation.m
% 频域实现二维相关（图案识别）：
%   f(x,y) (+) h(x,y) = IFFT{ F*(kx,ky) H(kx,ky) }      [书式 (1.2-7)(1.2-8)]
% 原书用两幅相同的 smiley.bmp（自相关，图 1.14 中心亮斑 = 匹配成功）；
% 把 I2 换成别的图案即复现图 1.15 的互相关（中心无亮斑）。
% smiley.bmp 可用画图板自制：256x256 黑底白色笑脸。
clear

I1=imread('smiley.bmp','bmp'); %Input image 1 (reference image)   % 参考图像
I1=I1(:,:,1);
figure(1) %displaying input image 1
colormap(gray(255));
image(I1)
axis off

FI1=fft2(I1);
max1=max(FI1);
max2=max(max1);
scale=1.0/max2;
FI1=FI1.*scale;                                                   % 归一化

I2=imread('smiley.bmp','bmp'); %Input image 2 (image to be recognized)
I2=I2(:,:,1);                                                     % 待识别图像
figure(2) %displaying input image 2
colormap(gray(255));
image(I2)
axis off

FI2=fft2(I2);
max1=max(FI2);
max2=max(max1);
scale=1.0/max2;
FI2=FI2.*scale;

% 说明：共轭放在哪个谱上只决定输出是否镜像翻转；实图像取 abs 显示后无差别。
FPR=FI1.*conj(FI2);%calculating correlation             % 频域相乘算相关
PR=ifft2(FPR);
PR=fftshift(PR);                                        % 相关峰搬到图像中心

max1=max(PR);
max2=max(max1);
scale=1.0/max2;
PR=PR.*scale;

figure(3)%display of correlation in spatial domain      % 相关结果灰度显示
colormap(gray(255));
image(abs(256*PR));
axis off
