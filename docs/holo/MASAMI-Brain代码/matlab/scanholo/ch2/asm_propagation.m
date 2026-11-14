% ===== 书 2.3 · 式 (2.3-5) 角谱传递函数 H(kx,ky;z) 的严格实现（ASM 法） =====
% 改编自你已有的 ASM.m（"严格对齐书本公式版"，规范化重写）：
% 一个振幅-相位物体（过冷水滴）经 z=5cm 传播生成同轴全息图，
% 再用 H(kx,ky;-z) 反向传播重建，验证角谱法"记录-重建"闭环。
%
% 核心公式（书中约定，前向变换为 +j 指数）：
%     H(kx,ky;z) = exp[ -j k0 z sqrt(1 - kx^2/k0^2 - ky^2/k0^2) ]   (2.3-5)
%     Psi(kx,ky;z) = Psi_p0(kx,ky) * H(kx,ky;z)                     (2.3-14)
%     psi_p(x,y;z) = IFT{ Psi(kx,ky;z) }                            (2.3-8)
%
% 差异提示：
% 1) 全书最大的坑：书的前向变换是 +j 指数，而 MATLAB 的 fft2 是 -j 指数。
%    严格按书公式实现前向变换必须写 ifft2(U)*N^2；逆变换写 fft2(A)/N^2。
%    若直接用 fft2 当"前向"，H 的指数符号要整体反号，否则传播方向颠倒
%    （强度图恰好看不出错，复振幅/相位全错）。
% 2) ifftshift 先把"中心在数组中央"的场搬到 DFT 原点排布，变换后 fftshift 搬回；
%    不做这一步相位会带 (-1)^(i+j) 锯齿。
% 3) 频率格间距用 DFT 真实值 df = 1/(N*dx)（不是书中 (N-1) 标定习惯，见 ch1 坑 3）。

clear
N = 512;                % 矩阵大小
lambda0 = 532e-9;       % 真空波长 lambda_0（m）
dx = 5e-6;              % 像素尺寸（m）
L = N*dx;               % 靶面尺寸（m）
z_prop = 0.05;          % 传播距离 z = 5 cm
k0 = 2*pi/lambda0;      % 波数 k_0

x = (-N/2:N/2-1)*dx;
y = (-N/2:N/2-1)*dx;
[X, Y] = meshgrid(x, y);

% ---- 虚拟构建"过冷水滴"物体（振幅 0.95 + 相位 1.2 rad 的圆盘） ----
U_obj = ones(N, N);
droplet_radius = 50e-6;
mask = (X.^2 + Y.^2) <= droplet_radius^2;
U_obj(mask) = 0.95*exp(1j*1.2);

% ---- 书中傅里叶变换协议封装（前向 +j = ifft2*N^2，逆向 -j = fft2/N^2） ----
% （函数体在文件末尾：ft2_book / ift2_book）

% ---- 频域角频率网格 (kx, ky) [rad/m] 与传递函数 (2.3-5) ----
df = 1/L;                       % DFT 频率格间距 [cycle/m]
fx = (-N/2:N/2-1)*df;
fy = (-N/2:N/2-1)*df;
[Fx, Fy] = meshgrid(fx, fy);
Kx = 2*pi*Fx;                   % 书中的 k_x
Ky = 2*pi*Fy;                   % 书中的 k_y

root_term = 1 - (Kx.^2 + Ky.^2)/k0^2;
root_term(root_term < 0) = 0;   % 滤除消逝波（k^2 > k0^2 的倏逝分量）

H_forward  = exp(-1j*k0*z_prop*sqrt(root_term));    % 式 (2.3-5)，z>0
H_backward = exp(-1j*k0*(-z_prop)*sqrt(root_term)); % 反向传播：z 取负

% ---- 正向传播：生成同轴（数字）全息图 ----
A_obj_0 = ft2_book(U_obj, N);           % 初始角谱 Psi_p0(kx,ky)
A_det   = A_obj_0 .* H_forward;         % 频域乘传递函数（式 2.3-14）
U_det   = ift2_book(A_det, N);          % 逆变换回空域（式 2.3-8）
I_hologram = abs(U_det).^2;             % 相机记录的强度（全息图）

% ---- 数字重构：反向传播 -z（取 sqrt(I) 当复振幅估计） ----
U_recon0 = sqrt(I_hologram);
A_recon  = ft2_book(U_recon0, N).*H_backward;
U_final  = ift2_book(A_recon, N);

% ---- 可视化 ----
figure(1)
subplot(2,2,1); imagesc(x*1e3, y*1e3, abs(U_obj)); colormap gray; axis image
title('Object amplitude (droplet, 0.95)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,2,2); imagesc(x*1e3, y*1e3, angle(U_obj)); colormap jet; axis image
title('Object phase (1.2 rad disc)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,2,3); imagesc(x*1e3, y*1e3, I_hologram); colormap gray; axis image
title('Inline hologram |U(z=5cm)|^2 (Fresnel rings)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,2,4); imagesc(x*1e3, y*1e3, abs(U_final)); colormap gray; axis image
title('Reconstructed amplitude (back-propagation)'); xlabel('x (mm)'); ylabel('y (mm)')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function A = ft2_book(U, N)
% 书式前向变换（+j 指数）：F(kx,ky) = 双重积分 U exp[j(kx x+ky y)] dx dy
A = fftshift(ifft2(ifftshift(U)))*N^2;
end

function U = ift2_book(A, N)
% 书式逆变换（-j 指数）：U(x,y) = (1/4pi^2) 双重积分 F exp[-j(kx x+ky y)] dkx dky
U = fftshift(fft2(ifftshift(A)))/N^2;
end
