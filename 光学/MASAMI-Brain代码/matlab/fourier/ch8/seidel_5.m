function[w]=seidel_5(u0,v0,X,Y,...
            wd,w040,w131,w222,w220,w311)
% ===== 书 8.2.2 节:seidel_5 前 5 项 Seidel 像差 + 离焦的波前 OPD =====
% 书 8.2.2 节列出的函数(原书第 1-30 行),实现式 (8.2):
%   W(u0; x,y) = wd*(x^2+y^2) + w040*(x^2+y^2)^2 + w131*u0*(x^2+y^2)*x
%              + w222*u0^2*x^2 + w220*u0^2*(x^2+y^2) + w311*u0^3*x
% 输入:
%   u0,v0 - 归一化像面坐标(标量;像点方向用 beta = atan2(v0,u0))
%   X,Y   - 归一化出瞳坐标网格(meshgrid 生成;物理坐标除以出瞳半径 wXP)
%   wd    - 离焦;w040 - 球差;w131 - 彗差;
%   w222  - 像散;w220 - 场曲;w311 - 畸变(单位:米,习惯用"波长数"x lambda)
% 要点:
%   - 旋转对称系统只需在子午面内写多项式:第 19-20 行把瞳坐标旋转 beta,
%     使 x' 轴指向像点方向(书 8.2.1 节);
%   - 五个初级像差满足 k+l = 4(下标 = 余弦幂/径向幂/场幂,式 8.1);
%   - 网格半径 >1 的"多余"值照常计算,后续用瞳函数乘 circ 时自然清零
%     (书原话:These extraneous values can be removed later)。
% MATLAB/Python 差异:Python 版在 python/fourier/ch8/ch8_aberrations.py
%   的 seidel_5();逐行对应,无索引差异(全为逐元素运算)。
% ------------- 以下为书原代码(第 1-30 行) -------------
% seidel_5
% Compute wavefront OPD for first 5 Seidel wavefront
% aberration coefficients + defocus
%
% u0,v0 - normalized image plane coordinate
% X,Y - normalized pupil coordinate arrays
%         (like from meshgrid)
% wd-defocus; w040-spherical; w131-coma;
% w222-astigmatism; w220-field curvature;
% w311-distortion

beta=atan2(v0,u0);     % image rotation angle
u0r=sqrt(u0^2+v0^2);   % image height

% rotate grid
Xr=X*cos(beta)+Y*sin(beta);
Yr=-X*sin(beta)+Y*cos(beta);

% Seidel polynomials
rho2=Xr.^2+Yr.^2;
w=wd*rho2+...
    w040*rho2.^2+...
    w131*u0r*rho2.*Xr+...
    w222*u0r^2*Xr.^2+...
    w220*u0r^2*rho2+...
    w311*u0r^3*Xr;
end
