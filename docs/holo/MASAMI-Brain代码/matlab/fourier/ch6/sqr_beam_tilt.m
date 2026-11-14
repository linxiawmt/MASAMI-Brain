% ===== 书 6.1 节:sqr_beam_tilt 倾斜方形光束例子(逐行忠实 + 中文注释) =====
% 书 6.1 节的做法:回到 5.3 节的 sqr_beam,Save As 成 sqr_beam_tilt,
% 在传播调用之前插入 4 行 tilt 代码(见下面"新增"块)。
% 参数与 5.3 节相同:L1=0.5 m、M=250、dx1=2 mm、w=0.051 m、
%   lambda=0.5 um、z=2000 m;新增 alpha=5.0e-5 rad、theta=45 deg。
% 图 6.3:观察面光斑沿 45 度方向平移 r = z*tan(alpha) = 0.1 m。
% 采样限制试验(书 6.1 节末):
%   - alpha = 7.6e-5(式 6.6 的判据值,B1 取 98 cycles/m):光束贴近阵列
%     边缘,幅度剖面出现 wrap-around(周期延拓把能量送进对面一侧);
%   - alpha = lambda/(2*dx1) = 1.25e-4:图案被网格边界切成两半。
% MATLAB/Python 差异:Python 版 ch6_sqr_beam_tilt.py;y=0 剖面在
%   I2(M/2+1,:),numpy 为 I2[M//2,:];rect 非基础 MATLAB 内置
%   (附录 C.3 有书自己的实现),无工具箱时可用逻辑掩码等价。
% ------------- 以下为书脚本(插入块按书 6.1 节) -------------
% sqr_beam_tilt tilted beam example
%
L1=0.5;               %side length
M=250;                %number of samples
dx1=L1/M;             %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;

lambda=0.5*10^-6;     %wavelength
k=2*pi/lambda;        %wavenumber
w=0.051;              %source half width (m)
z=2000;               %propagation dist (m)

[X1,Y1]=meshgrid(x1,y1);
u1=rect(X1/(2*w)).*rect(Y1/(2*w)); %src field
I1=abs(u1.^2);        %src irradiance
%
figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m');

% ===== 书 6.1 节新增:传播调用前插入 tilt =====
deg=pi/180;
alpha=5.0e-5; %rad
theta=45*deg;
[u1]=tilt(u1,L1,lambda,alpha,theta);

u2=propTF(u1,L1,lambda,z); %propagation

x2=x1;                %obs coords
y2=y1;
I2=abs(u2.^2);        %obs irrad

figure(2)             %display obs irrad
imagesc(x2,y2,I2);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);
%
figure(3)             %irradiance profile
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
title(['z= ',num2str(z),' m']);
%
figure(4)             %plot obs field mag
plot(x2,abs(u2(M/2+1,:)));
xlabel('x (m)'); ylabel('Magnitude');
title(['z= ',num2str(z),' m']);
%
figure(5)             %plot obs field phase
plot(x2,unwrap(angle(u2(M/2+1,:))));
xlabel('x (m)'); ylabel('Phase (rad)');
title(['z= ',num2str(z),' m']);
