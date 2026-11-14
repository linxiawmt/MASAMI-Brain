% ===== 书 5.5 节:sqr_beam_FF 夫琅禾费例子(书"Save as"改动版) =====
% 书 5.5 节:把 sqr_beam 用"Save as"存成 sqr_beam_FF 后做如下改动:
%   w=0.011;              %source half width (m)
%   [u2,L2]=propFF(u1,L1,lambda,z);
%   dx2=L2/M;
%   x2=-L2/2:dx2:L2/2-dx2; %obs coords
%   y2=x2;
%   imagesc(x2,y2,nthroot(I2,3)); %stretch image contrast
% 参数:w = 0.011 m、z = 2000 m -> NF = w^2/(lambda z) = 0.12
% (书:"reasonable for the Fraunhofer approximation")。
% L2 = lambda*z/dx1 = 0.5 m、dx2 = lambda*z/L1 = 2 mm(本例恰为临界采样)。
% 注意:书改动清单里 imagesc 直接用了 I2——I2 = abs(u2.^2) 必须先算
% (原书没列这一行,补在第 26 行位置)。
% 依赖:同文件夹 propFF.m(书 5.5 节)。
% MATLAB/Python 差异:Python 版 ch5_fraunhofer.py;nthroot(I2,3) 对应
%   I2**(1/3)。
% ------------- 以下为改动后的完整脚本 -------------
L1=0.5;               %side length
M=250;                %number of samples
dx1=L1/M;             %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;

lambda=0.5*10^-6;     %wavelength
k=2*pi/lambda;        %wavenumber
w=0.011;              %source half width (m)   <- 书改动行
z=2000;               %propagation dist (m)

[X1,Y1]=meshgrid(x1,y1);
u1=rect(X1/(2*w)).*rect(Y1/(2*w)); %src field
I1=abs(u1.^2);        %src irradiance

[u2,L2]=propFF(u1,L1,lambda,z);   % <- 书改动行(Fraunhofer)

dx2=L2/M;             % <- 书改动行
x2=-L2/2:dx2:L2/2-dx2; %obs coords  <- 书改动行
y2=x2;
I2=abs(u2.^2);        %obs irrad(原书改动清单省略了这行,补上)

figure(1)             %irradiance image(1/3 次根拉出旁瓣)
imagesc(x2,y2,nthroot(I2,3)); %stretch image contrast <- 书改动行
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);

figure(2)             %irradiance profile
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
title(['z= ',num2str(z),' m']);

% ---------- 数值检查(对照式 (5.22) 解析) ----------
peak_an = (4*w^2/(lambda*z))^2;       % 式 (5.22) 峰值
null_an = lambda*z/(2*w);             % sinc 第一零点:2w x/(lambda z) = 1
side_an = sinc(1.5)^2 * peak_an;      % 第一旁瓣(sinc 为信号处理工具箱定义)
fprintf('NF = w^2/(lambda z) = %.3f(书:0.12)\n', w^2/(lambda*z));
fprintf('L2 = %.3f m,dx2 = %.1f mm(临界采样)\n', L2, dx2*1e3);
fprintf('I2 峰值:解析 %.4f,实测 %.4f\n', peak_an, max(I2(:)));
fprintf('sinc 第一零点 x = lambda z/(2w) = %.4f m\n', null_an);
fprintf('第一旁瓣:解析 %.4f(0.045*I2peak),实测 %.4f\n', ...
        side_an, max(I2(M/2+1, M/2+40:M/2+90)));
fprintf('总功率:源面 %.6f -> 观察面 %.6f\n', ...
        sum(sum(I1))*dx1^2, sum(sum(I2))*dx2^2);
