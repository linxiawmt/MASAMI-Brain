% ===== 拓展实战（编者，原书无）· 同轴数字全息的背景相除法 =====
% 改编自你已有的 holography.m / Twins_Image.m（规范化重写）：
% 过冷水滴（振幅 0.95 + 相位 1.2 rad 的圆盘）经 z=5cm 角谱传播生成同轴全息图，
% 相机面叠加高斯激光背景与两粒灰尘 -> 背景相除 I_raw/I_bg -> sqrt 估计复振幅 -> 反向传播重建。
%
% 与本章 OSH 的关系（笔记 3.5 节"孪生像问题"对照）：
% - 背景相除只能去掉乘性背景（灰尘/照明不均），对同轴全息的孪生像虚影无能为力——
%   重建振幅图里水滴周围仍有发散的"孪生像波纹"。
% - OSH 路线（本章）：一次 2-D 扫描同时记录 sin/cos 两张全息图，复数相加 H_cos + j*H_sin
%   即可在同轴记录下彻底消除孪生像（式 3.5-5）。
%
% 差异提示：
% 1) 全书最大的坑：书的前向变换是 +j 指数，严格实现用 ifft2*N^2；逆变换用 fft2/N^2（文件末尾封装）。
%    你旧版 holography.m 的 ifft2(...)*N^2 + fft2(...)/N^2 组合正是书协议（与 2.3 节 asm_propagation 一致）。
% 2) 连通域定位用 regionprops（Image Processing Toolbox）；去 toolbox 的简化版见 Python 版。

clear
N = 512; lambda0 = 532e-9; dx = 5e-6; L = N*dx; z_prop = 0.05; k0 = 2*pi/lambda0;
x = (-N/2:N/2-1)*dx; y = (-N/2:N/2-1)*dx;
[XX, YY] = meshgrid(x, y);

% ---- 1. 模拟背景：高斯激光 + 两粒灰尘 ----
w0 = 1.5e-3;
I_laser = exp(-2*(XX.^2+YY.^2)/w0^2);
t_noise = ones(N);
dust1 = ((XX-0.5e-3).^2 + (YY-0.5e-3).^2) <= (20e-6)^2;
dust2 = ((XX+0.6e-3).^2 + (YY-0.4e-3).^2) <= (15e-6)^2;
t_noise(dust1 | dust2) = 0.2;
I_bg = I_laser.*t_noise;

% ---- 2. 原始物体：过冷水滴（振幅 0.95 + 相位 1.2 rad） ----
U_obj = ones(N);
mask = XX.^2 + YY.^2 <= (50e-6)^2;                 % 半径 50 um
U_obj(mask) = 0.95*exp(1j*1.2);

% ---- 3. 角谱正/反向传播（书式协议，自动滤消逝波） ----
df = 1/L; fx = (-N/2:N/2-1)*df; fy = (-N/2:N/2-1)*df;
[KY, KX] = meshgrid(fy, fx);
root = 1 - (KX.^2+KY.^2)/k0^2; root(root < 0) = 0;   % 滤除消逝波
H_fwd = exp(-1j*k0*z_prop*sqrt(root));
H_bwd = exp(+1j*k0*z_prop*sqrt(root));               % 反向传播：z 取 -z_prop

U_det = ift2_book(ft2_book(U_obj).*H_fwd);           % 传播到相机面
I_raw = (abs(U_det).^2).*I_bg;                       % 相机记录（含背景）
I_proc = I_raw./I_bg;                                % 背景相除（乘性噪声消除）
U_final = ift2_book(ft2_book(sqrt(I_proc)).*H_bwd);  % sqrt 强度当复振幅估计，反向传播

% ---- 4. 定位分析 ----
amp = abs(U_final);
amp_n = (amp - min(amp(:)))/(max(amp(:)) - min(amp(:)));
binary = amp_n < 0.45;
stats = regionprops(binary, 'Centroid', 'Area');
fprintf('检测到 %d 个连通域；', length(stats));
for kk = 1:length(stats)
    c = stats(kk).Centroid;
    d_eq = 2*dx*sqrt(stats(kk).Area/pi)*1e6;
    fprintf('质心 (%+.3f, %+.3f) mm, 等效直径 %.1f um;  ', ...
            (c(1)-N/2)*dx*1e3, (c(2)-N/2)*dx*1e3, d_eq);
end
fprintf('\n注意：重建振幅中水滴周围仍可见孪生像虚影波纹——背景相除无法消除孪生像。\n');

% ---- 5. 可视化 ----
figure('Position', [100, 100, 1200, 800]);
subplot(2,3,1); imagesc(x*1e3, y*1e3, abs(U_obj)); colormap gray; axis image; colorbar
title('1. Object amplitude |U_obj| (z=0)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,3,2); imagesc(x*1e3, y*1e3, I_raw); axis image; colorbar
title('2. Raw hologram (laser bg + dust)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,3,3); imagesc(x*1e3, y*1e3, I_proc); axis image; colorbar
title('3. Background-divided hologram'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,3,4); imagesc(x*1e3, y*1e3, amp); axis image; colorbar
title('4. Reconstructed amplitude (twin-image haze)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,3,5); imagesc(x*1e3, y*1e3, angle(U_final)); colormap jet; axis image; colorbar
title('5. Reconstructed phase (rad)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(2,3,6); imagesc(x*1e3, y*1e3, amp); axis image; hold on
title('6. Detection'); xlabel('x (mm)'); ylabel('y (mm)')
for kk = 1:length(stats)
    c = stats(kk).Centroid;
    plot((c(1)-N/2)*dx*1e3, (c(2)-N/2)*dx*1e3, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
hold off

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function F = ft2_book(U)
% 书式前向变换（+j 指数）
N = size(U, 1);
F = fftshift(ifft2(ifftshift(U)))*N^2;
end

function U = ift2_book(A)
% 书式逆变换（-j 指数）
N = size(A, 1);
U = fftshift(fft2(ifftshift(A)))/N^2;
end
