% ===== 书 5.3 节:sqr_beam 方形光束传播例子(逐行忠实 + 中文注释) =====
% 书 5.3 节 New M-file(命名为 "sqr_beam")的主脚本(原书第 1-47 行)。
% 复现图 5.1(源面)、图 5.2(观察面辐照度)、图 5.3(幅度/相位剖面)。
% 参数:L1 = 0.5 m、M = 250、dx1 = 2 mm、w = 0.051 m(跨孔径 51 样本)、
%   lambda = 0.5 um、z = 2000 m、NF = w^2/(lambda z) = 2.6(书:"reasonable")。
% 依赖:同文件夹 propTF.m(书 5.1 节)/ propIR.m(书 5.2 节)。
% MATLAB/Python 差异:Python 版 ch5_sqr_beam.py;imagesc + axis xy 对应
%   imshow(origin='lower');y=0 剖面在 I2(M/2+1,:),numpy 在 I2[M//2,:]。
% ------------- 以下为书原代码(第 1-47 行,22 行的 propIR 变体见文末) -------------
% sqr_beam propagation example
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

% ---------- 数值检查(对照书 5.3 节的口算与图 5.2/5.3) ----------
% 注意:rect 在基础 MATLAB 里不是内置函数(信号处理工具箱才有),
% 若无工具箱,书第 14 行可用逻辑掩码等价实现:
%   u1 = (abs(X1) < w) .* (abs(Y1) < w);
fprintf('dx1 = %.1f mm(书:2 mm)\n', dx1*1e3);
fprintf('跨孔径样本数 = 2w/dx1 = %.0f(书:51)\n', 2*w/dx1);
fprintf('NF = w^2/(lambda z) = %.2f(书:2.6)\n', w^2/(lambda*z));
I2_tf = I2;
u2_ir = propIR(u1,L1,lambda,z);   % 书 5.3 节末尾:把 22 行换成 propIR 再跑
fprintf('max|I2_TF - I2_IR|/max(I2) = %.2e(临界采样下应几乎为 0)\n', ...
        max(max(abs(I2_tf - abs(u2_ir.^2))))/max(I2_tf(:)));
fprintf('I2 峰值 = %.4f(图 5.2(b) 约 1.95;x=0 中央峰 1.92,±0.01 谷 0.89)\n', ...
        max(I2_tf(:)));
fprintf('总功率:源面 %.6f -> 观察面 %.6f(习题 5.3 功率守恒)\n', ...
        sum(sum(I1))*dx1^2, sum(sum(I2_tf))*dx1^2);
