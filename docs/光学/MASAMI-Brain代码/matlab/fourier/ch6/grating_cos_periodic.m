% ===== 书 6.4.4 节:grating_cos 无限周期模型(式 6.28 + 分数功率) =====
% 把光栅结构铺满整个阵列(D1 = L1),且整数个周期正好跨满阵列:
%   L1/P = integer(式 6.28)。这里 D1=L1=1e-2,L1/P=100(100 个周期)。
% DFT 的周期延拓让"阵列边界完美连续"的源等价于无限周期结构——
% 夫琅禾费图样退化为衍射阶位置上的样本 delta 函数。
% 分数功率(书给的代码,假定单位源场幅度):
%   Fpwr = I2*dx2*dx2/(L1^2);
% 结果(图 6.16):零阶 25%、±1 一阶各 6.25%。
% 注意(笔记勘误):这组数字是**入射功率**的分数,不是书正文说的
%   "除以光栅透射的总功率"——振幅光栅 1/2*(1-cos) 的透射率平方均值
%   只有 3/8,若除以真实透射功率,零阶应占 2/3、±1 阶各 1/6。
%   书代码的分母 L1^2 = 单位入射场在光栅面积上的功率,与正文措辞矛盾。
% MATLAB/Python 差异:Python 版 ch6_gratings.py 的 periodic model 部分。
% ------------- 以下为本笔记脚本(书 6.4.4 节改法) -------------
% grating_cos_periodic infinite periodic model
lambda=0.5e-6;  %wavelength
f=0.5;          %focal distance
P=1e-4;         %grating period

L1=1e-2;        %array side length = grating side length(式 6.28)
D1=L1;          %光栅铺满阵列
M=500;          %# samples
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1;
[X1,Y1]=meshgrid(x1,x1);

u1=1/2*(1-cos(2*pi*X1/P)).*rect(X1/D1).*rect(Y1/D1);

[u2,L2]=propFF(u1,L1,lambda,f);
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; y2=x2;
I2=abs(u2).^2;

Fpwr=I2*dx2*dx2/(L1^2);      %书原代码:每个样本的功率分数
%
figure(1)                    %图 6.16:中心行切片
plot(x2,Fpwr(M/2+1,:));
xlabel('x (m)'); ylabel('Fractional Power');
% 阶位置:零阶 x=0 处 Fpwr=0.25;±lambda*f/P=±2.5 mm 处各 0.0625
