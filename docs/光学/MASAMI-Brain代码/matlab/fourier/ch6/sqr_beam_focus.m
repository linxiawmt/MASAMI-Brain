% ===== 书 6.2 节:sqr_beam_focus 聚焦方形光束例子(逐行忠实 + 中文注释) =====
% 书 6.2 节的做法:回到 sqr_beam,Save As 成 sqr_beam_focus,在传播调用
% 之前插入 focus 代码(见"新增"块)。
% 参数与 5.3 节相同;zf = 2000 m 恰等于传播距离 z,观察面上是一个小光斑
% ——实际上就是一幅缩放了的夫琅禾费图样(书 6.2 节原话,见 6.3 节)。
% 边界试验(书 6.2 节末,B1 = 98 cycles/m 时式 6.11 给 |zf| >= 671 m):
%   - zf = -671:发散波,观察面上光束明显推过阵列边界;
%   - zf = +671:光束先聚焦再散开,图样变小、看似合理,但阵列边缘处
%     相位混叠误差刚刚开始出现(看 unwrap 相位剖面)。
% MATLAB/Python 差异:Python 版 ch6_sqr_beam_focus.py;propFF 对比实验
%   (聚焦图样 = 缩放夫琅禾费)见 Python 脚本与笔记正文。
% ------------- 以下为书脚本(插入块按书 6.2 节) -------------
% sqr_beam_focus focused beam example
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

% ===== 书 6.2 节新增:传播调用前插入 focus =====
zf=2000;
[u1]=focus(u1,L1,lambda,zf);

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

% ===== 书 6.2 节末的边界试验(把 zf 换成下列值重跑) =====
% zf=-671;   % 发散:光束推过观察阵列边界
% zf=+671;   % 先聚焦:图样变小,边缘刚开始出现相位混叠
