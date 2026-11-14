% ===== 书 1.1 例 1.1 · Table 1.3 fft2Drect.m =====
%fft2Drect.m %Simulation of Fourier transformation of a 2-D rect function
% 二维 rect 函数（M x M 正方形）傅里叶变换仿真，验证
%   rect(x/a, y/b) <=> ab*sinc(a*kx/2pi, b*ky/2pi)
% 运行后按提示输入奇数 M（3~33）；书中示例 M=11：
%   a = 11/256 = 0.0429，第一零点 kx0 = 2*pi/a = 146.23 rad/(unit of length)
clear
L=1; %display area is L by L, L has unit of length   % 显示区域 L x L
N=256; % number of sampling points                   % 采样点数
dx=L/(N-1); % dx : step size                         % 采样步长

% Create square image, M by M square, rect(x/a,y/a), M=odd number
M=input('M (size of rect(x/a,y/a), enter odd numbers from 3-33-)=');
a=M/256;
kx0=2*pi/a;
sprintf('a= %0.5g[unit of length]',a)
sprintf('kx0 (first zero)= %0.5g[radian/unit of length]',kx0)
R=zeros(256); %assign a matrix (256x256) of zeros    % 背景全 0
r=ones(M); % assign a matrix (MxM) of ones           % 中间全 1 的 M x M 方块
n=(M-1)/2;
R(128-n:128+n,128-n:128+n)=r;                        % 把方块嵌进 256 网格正中
%End of creating square input image M by M

%Axis Scaling                                        % 坐标轴标定（零频移到中央）
for k=1:256
    X(k)=1/255*(k-1)-L/2;                            % 空间坐标，范围 [-0.5, 0.5]
    Y(k)=1/255*(k-1)-L/2;
     %Kx=(2*pi*k)/((N-1)*dx)
    %in our case, N=256, dx=1/255
    Kx(k)=(2*pi*(k-1))/((N-1)*dx)-((2*pi*(256-1))/((N-1)*dx))/2;  % 角频率轴
    Ky(k)=(2*pi*(k-1))/((N-1)*dx)-((2*pi*(256-1))/((N-1)*dx))/2;
end

%Image of the rect function                          % 图 1：rect 灰度图
figure(1)
image(X+dx/2,Y+dx/2,255*R);
title('rect function: gray-scale plot')
xlabel('x')
ylabel('y')
colormap(gray(256));
axis square

%Computing Fourier transform                         % 二维 FFT + fftshift 中心化
% 说明：书的前向变换约定是"正指数"，严格写法是 ifft2(R)*256^2；
% rect 是实偶函数，fft2(R)/256^2 与之完全相等，故书里直接用 fft2。
FR=(1/256)^2*fft2(R);
FR=fftshift(FR);

% plot of cross-section of rect function             % 图 2：rect 的中心横截面
figure(2)
plot(X+dx/2,R(:,127))
title('rect function: cross-section plot')
xlabel('x')
ylabel('rect(x/a)')
grid
axis([-0.5 0.5 -0.1 1.2])

%Centering the axis and plot of cross-section of transform along kx
% 图 3：|F| 沿 kx 的横截面（幅度 x10 显示；横轴再平移半个频格对准零频）
figure(3)
plot(Kx-pi/(dx*(N-1)),10*abs(FR(:,127)))
title('Square-absolute value of Fourier transform of rect function: cross-section plot')
xlabel('kx')
ylabel('|a*b*sinc(a*kx/2pi)|')
axis([-800 800 0 max(max(abs(FR)))*10.1])
grid

%Mesh the Fourier transformation                     % 图 4：|F|^2 三维网格图
figure(4);
mesh(Kx,Ky,(abs(FR)).^2)
title('Square-absolute value of Fourier transform of rect function: 3-D plot,scale arbitrary')
xlabel('kx')
ylabel('ky')
axis square

%Image of the Fourier transformation of rectangular function
% 图 5：|F|^2 灰度图（gain=10000 压动态范围，暗背景上才看得清亮斑）
figure(5)
gain=10000;
image(Kx,Ky,gain*(abs(FR)).^2/max(max(abs(FR))).^2)
title('Square-absolute value of Fourier transform of the rect function: gray-scale plot')
xlabel('kx')
ylabel('ky')
axis square
colormap(gray(256))
