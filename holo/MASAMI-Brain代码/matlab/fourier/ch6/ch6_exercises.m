% ===== 书 6.5 节习题 6.2/6.3/6.5/6.6/6.7:MATLAB 骨架 =====
% 数值结果与图见笔记正文与 python/fourier/ch6/ch6_exercises.py 实跑。
% 依赖:propTF.m(ch5)、tilt.m、focus.m、ucomb.m(本目录)。

%% ---- 习题 6.2:tilt + focus 同时作用(sqr_beam,TF,z=2000 m) ----
% 参数:alpha=5.0e-5 rad、theta=45 deg、zf=4000 m。
% 做法:在传播调用前依次乘两个透过率(乘法可交换):
%   [u1]=tilt(u1,L1,lambda,alpha,theta);
%   [u1]=focus(u1,L1,lambda,zf);
%   u2=propTF(u1,L1,lambda,z);
% 预期:光斑沿 45 度方向平移 z*tan(alpha)=0.1 m,且因 zf>z(先聚焦后
% 观察而未到焦点……注意 zf=4000 > z=2000,观察面在焦点之前),斑比
% 纯聚焦(zf=z)大、比无聚焦小。

%% ---- 习题 6.3:柱面透镜聚焦(只在一个横向轴聚焦) ----
% (a) 把 focus 的 chirp 只作用在 x 上(实现见本文件末尾的 focus_cyl):
%     [u1]=focus_cyl(u1,L1,lambda,2000);
% (b) 用 sqr_beam 演示:zf=2000 m——观察面上 x 方向收到小斑、
%     y 方向仍是衍射展宽(一个"一条线"的焦线)。
% (c) 柱面透镜焦面的场只有一个轴是傅里叶变换(另一个轴仍带传播/菲涅耳
%     行为),不能像式 (6.18)/(6.19) 那样整体写成二维夫琅禾费图样。

%% ---- 习题 6.5:波带片(cos 型与二值型) ----
% t_A(r) = 1/2*(1+cos(k*r^2/(2*f)))*circ(r/w)   (式 6.29)
% t_A(r) = 1/2*(1+sign(cos(k*r^2/(2*f))))*circ(r/w)   (式 6.30)
% 参数:L1=50 mm、M=500、lambda=0.5 um、w=6.25 mm、f=10 m;z=5/8/10/12 m。
% circ 约定:式 (6.29) 按附录 C.1(out=abs(r)<=1),半径 = w;
% 用逻辑掩码等价:r2=sqrt(X1.^2+Y1.^2); P=double(r2<=w);
% k=2*pi/lambda;
% tA1=0.5*(1+cos(0.5*k*r2.^2/f)).*P;          %cos 型(0.5*k/f = k/(2f))
% tA2=0.5*(1+sign(cos(0.5*k*r2.^2/f))).*P;    %二值型(MATLAB 自带 sign)
% for z=[5 8 10 12]
%     u2=propTF(tA1,L1,lambda,z);  %z=10 m 处聚焦(放大图样)
%     ...显示 I2=abs(u2).^2 与剖面
% end
% (g) z=10 m 焦点峰值:cos 型 vs 二值型——二值型把相位钳到 0/pi,
%     等效"全切全偏",一级效率更高(峰值对比见 Python 实跑);
% (h) 与同口径(f=10 m,半径 w)的理想正透镜焦斑比:透镜斑远小而高
%     ——波带片有高阶焦点与残余零阶,只在小 NA/特殊波段下可用。

%% ---- 习题 6.6:周期函数采样敏感性(grating_sqr 改 M) ----
% P/dx1 = M*P/L1 必须是整数(每周期整数个样本,ucomb 才能落位)。
% grating_sqr:P=1e-4、L1=1e-2 -> 每周期样本数 = M/100。
% M=1002:每周期 10.02 个样本,delta 落位漂移,夫琅禾费阶"发毛";
% M=1000 之上下一个可用值是 1100(每周期 11 个样本)。

%% ---- 习题 6.7:相位光栅(改 grating_cos / grating_sqr) ----
% (a) u1=exp(j*m*cos(2*pi*X1/P)).*rect(X1/D1).*rect(Y1/D1);
% (b) u1=exp(j*m*(ux)).*rect(X1/D1).*rect(Y1/D1);   %ux 同 grating_sqr(方波)
% (c) u1=exp(j*m*(ut)).*rect(X1/D1).*rect(Y1/D1);   %ut:三角波 rect 换 tri
% m 取 pi、pi/2、pi/4:纯相位光栅不吸收(总透射恒为 1),零阶功率
% J0(m)^2 随 m 减小、能量往高阶搬;J0 的第一个零点在 m≈2.405,
% 此时零阶完全消失(m=pi 时 J0(pi)^2=0.0926,还剩约 9%)。

%% ---- 本文件末尾:脚本局部函数(R2016b+ 支持) ----
function[uout]=focus_cyl(uin,L,lambda,zf)
% 习题 6.3(a):柱面聚焦——chirp 只作用在 x 轴,y 方向不变
[M,N]=size(uin);
dx=L/M; k=2*pi/lambda;
x=-L/2:dx:L/2-dx;
[X,~]=meshgrid(x,x);          %meshgrid 的第一输出:沿列(第 2 维)变化的 x 坐标
uout=uin.*exp(-j*k/(2*zf)*(X.^2));
end
