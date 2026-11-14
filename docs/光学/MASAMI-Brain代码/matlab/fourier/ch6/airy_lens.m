% ===== 书 6.3 节 / 习题 6.4:f/10 理想透镜焦面的 Airy 图样 =====
% 书 6.3 节:焦面场 = 入射场的夫琅禾费图样(z 换成 f,式 6.18),
% 用上一章的 propFF 直接算(z 参数填 f),不必用菲涅耳传播器硬扛
% 式 (6.15) 的采样要求。参数(习题 6.4(a)):f=250 mm、2wL=25 mm、
% M=250、L=250 mm、lambda=0.5 um;图 6.7(a) 用四次方根拉伸对比度。
% 解析对照(Airy 图样):I2(r) = (pi*wL^2/(lambda*f))^2 * [2*J1(v)/v]^2,
%   v = k*wL*r/f = pi*D*r/(lambda*f),D = 2wL;峰值 I2(0) = (pi*wL^2/(lambda*f))^2。
% 注意 circ 约定(式 6.13 按附录 C.1 / Goodman 约定,半径 = 分母):
%   P = circ(sqrt(x^2+y^2)/wL) 表示半径 wL 的圆孔。基础 MATLAB 没有
%   circ,这里用逻辑掩码等价实现(等价于 base.circle(X,Y,wL))。
% MATLAB/Python 差异:Python 版 ch6_lens_airy.py(用 base.circle)。
% ------------- 以下为本笔记脚本(对应书习题 6.4 / 图 6.7) -------------
% airy_lens focal-plane Airy pattern of an ideal f/10 lens
lambda=0.5e-6;        %wavelength
f=0.25;               %focal length (250 mm)
wL=0.0125;            %lens aperture radius (2wL=25 mm, f/10)
M=250;                %number of samples
L=0.25;               %array side length (250 mm)
dx=L/M;               %sample interval
x=-L/2:dx:L/2-dx;     %coords
[X,Y]=meshgrid(x,x);

u1=ones(M,M);                     %unit amplitude plane wave
P=double(sqrt(X.^2+Y.^2)<=wL);    %pupil: circ(sqrt(x^2+y^2)/wL), 半径 wL
u1=u1.*P;                         %lens pupil clips the plane wave

[u2,L2]=propFF(u1,L,lambda,f);    %Fraunhofer with z = f(式 6.18)
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; y2=x2;     %obs coords
I2=abs(u2).^2;
%
figure(1)                         %图 6.7(a):四次方根拉伸
imagesc(x2,y2,nthroot(I2,4));
axis square; axis xy; colormap('gray');
xlabel('x (m)'); ylabel('y (m)');
%
figure(2)                         %图 6.7(b):x 剖面
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
%
% analytic Airy pattern(习题 6.4(a) 的解析结果)
lf=lambda*f;
v=k*wL*x2/f;                      % = pi*D*x2/(lambda*f)
J1=besselj(1,v);
I2a=(pi*wL^2/lf)^2*(2*J1./v).^2;  %v=0 处 2*J1./v 的极限为 1
I2a(M/2+1)=(pi*wL^2/lf)^2;        %x2=0 样本(第 M/2+1 号)直接取峰值
figure(3)
plot(x2,I2a,'--');                %虚线解析曲线,与 figure(2) 实线对照
xlabel('x (m)'); ylabel('Irradiance');
