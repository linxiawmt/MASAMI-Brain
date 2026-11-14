function[uout]=tilt(uin,L,lambda,alpha,theta)
% ===== 书 6.1 节:tilt 倾斜透过率函数(逐行忠实 + 中文注释) =====
% 书 6.1 节列出的函数(原书第 1-20 行),实现式 (6.3):
%   t_a(x,y) = exp[jk(x cos(theta) + y sin(theta)) tan(alpha)]
% 即给波前加一个线性相位斜面,使光束朝 (cos(theta), sin(theta)) 方向
% 偏转 alpha 角(瞄准点 r = z tan alpha,图 6.2)。
% 要点:
%   - theta 从 x 轴量起(书注释 "x axis 0");式 (6.2) 中 theta = atan(y/x)
%     只是"径向方向"的定义,函数里 theta 是自由参数;
%   - 第 18-19 行用三个连续的句点(...)续行,一条表达式拆成两行;
%   - 倾斜角受采样限制(式 6.6):alpha <= lambda*(1/(2*dx) - B1),
%     倾斜本质是把源频谱平移 tan(alpha)/lambda(式 6.4/6.5)。
% MATLAB/Python 差异:Python 版在 python/fourier/ch6/ch6_transmittance.py
%   的 tilt();坐标轴生成见公共模块 base.coords(Python 侧)。
% ------------- 以下为书原代码(第 1-20 行) -------------
% tilt phasefront
% uniform sampling assumed
% uin - input field
% L - side length
% lambda - wavelength
% alpha - tilt angle
% theta - rotation angle (x axis 0)
% uout - output field

[M,N]=size(uin);        %get input field array size
dx=L/M;                 %sample interval
k=2*pi/lambda;          %wavenumber

x=-L/2:dx:L/2-dx;       %coords
[X,Y]=meshgrid(x,x);

uout=uin.*exp(j*k*(X*cos(theta)+Y*sin(theta))...
   *tan(alpha));       %apply tilt
end
