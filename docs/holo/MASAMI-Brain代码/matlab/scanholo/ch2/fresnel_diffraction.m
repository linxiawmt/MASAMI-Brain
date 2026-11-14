% ===== 书 2.3.2 · Table 2.1 Fresnel_diffraction.m（规范化重写） =====
% 正方形孔径的菲涅耳衍射仿真（卷积路线，式 (2.3-12)）：
%     psi_p(x,y;z) = psi_p0(x,y) * h(x,y;z)
% 输入：M x M 正方形孔 rect(x/a, y/a)，a = M/256 = 0.4336 cm
% 光源：He-Ne 红光 632.8 nm，观察屏 z = 5 cm（复现原书图 2.11）
%
% 与第 1 章拓展 1.2.3 的 fresnel_lsi.m 是同一物理的"书中正式版"。
%
% 差异提示：本脚本用"显式构造 h 再频域相乘"的卷积路线——卷积定理只依赖同一对
% 正/逆变换，与 +j/-j 指数约定无关，所以书中 fft2 写法的结果就是书的卷积
% （h 本身已按式 (2.3-11) 用 -j 二次相位构造）。"约定坑"只在用解析传递函数
% H(kx,ky;z)（式 2.3-13）时才咬人，见 asm_propagation.m。
% exp(-4*200*r^2) 是高斯窗：把脉冲响应截断在显示窗口内，抑制卷绕误差（原书写法）。

clear
L = 1;                  % L : length of display area（显示区域 L x L，单位 cm）
N = 256;                % N : number of sampling points（采样点数）
dx = L/(N-1);           % dx : step size（采样步长 1/255 cm）
M = 111;                % 正方形孔的像素宽度（a = M/256 = 0.4336 cm）
lambda0 = 0.6328e-4;    % He-Ne 红光波长（cm）
z = 5;                  % 观察距离（cm）

% ---- 输入：中心 M x M 正方形孔 rect(x/a, y/a) ----
R = zeros(256);                 % assign a matrix (256x256) of zeros
r = ones(M);                    % assign a matrix (MxM) of ones
n = (M-1)/2;
R(128-n:128+n, 128-n:128+n) = r;% 背景置 0、中心置 1 的方孔

% ---- 坐标轴标定（与书中 fft2Drect.m 相同的方式） ----
for k = 1:256
    X(k) = 1/255*(k-1) - L/2;   % 空间坐标 [-0.5, 0.5] cm
    Y(k) = 1/255*(k-1) - L/2;
end

% ---- 菲涅耳脉冲响应 h(x,y;z)（式 2.3-11 + 高斯窗） ----
sigma = pi/(lambda0*z);         % sigma = pi/(lambda*z)，二次相位系数
for rr = 1:256
    for cc = 1:256
        % j*(sigma/pi) = j/(lambda*z) 振幅因子（全局相位 exp(-j k0 z) 未算）
        h(rr,cc) = j*(sigma/pi)*exp(-4*200*(X(rr)^2+Y(cc)^2)) ...
                                *exp(-j*sigma*(X(rr)^2+Y(cc)^2));
    end
end

% ---- 频域实现卷积（照抄原书流程：FR.*H 后 ifft2 再 fftshift） ----
FR = (1/256)^2*fft2(R);         % 原书：FR=(1/256)^2*fft2(R)
H  = (1/256)^2*fft2(h);         % 原书：H=(1/256)^2*fft2(h)
hr = ifft2(FR.*H);
hr = (256^2)*hr;                % 原书：hr=(256^2)*hr
hr = fftshift(hr);              % 把卷积中心搬回图像中央
I  = abs(hr).^2;                % 光强

% ---- 可视化（对应原书 figure(1)(2)(3) 与 CCD 视图） ----
figure(1)
image(X, Y, 255*R); colormap(gray(256)); axis square; axis image
xlabel('cm'); ylabel('cm')
title('Input square aperture rect(x/a, y/a)')

figure(2)
plot(X+dx/2, R(:,127), '-o', 'MarkerSize', 2)   % 方孔中心横截面
grid; axis([-0.5 0.5 -0.1 1.2])
xlabel('cm'); title('Square aperture: cross-section')

figure(3)
plot(X+dx/2, abs(hr(:,127)))    % 衍射幅度中心横截面（原图 2.11c）
grid; xlabel('cm'); title('|psi_p(x,y;z=5cm)| cross-section')

figure(4)
imagesc(X, Y, I); axis image; colormap(gray(256)); colorbar
xlabel('x (cm)'); ylabel('y (cm)')
title('2-D diffraction pattern (CCD view) at z = 5 cm')
