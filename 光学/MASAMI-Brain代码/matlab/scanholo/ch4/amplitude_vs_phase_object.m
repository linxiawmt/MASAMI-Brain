% ===== 拓展实战（编者，原书无）· 同轴数字全息下 振幅物体 vs 相位物体 的成像对比 =====
% 改编自你已有的 Amplitude_and_Phase_Object.m（规范化重写）：
%   左：纯振幅粒子（强吸收 0.1，无相位）；右：纯相位粒子（无吸收，相位延迟 1.5 rad）
%   -> 角谱正向传播 5 cm 生成同轴全息图 -> 相机面取强度 -> sqrt 强度当复振幅估计 -> 反向传播重建。
%
% 与 4.1 节的关系（为什么扫描全息显微重要）：
% - 生物样品多为"相位物体"：几乎不吸收，只推迟相位——普通强度成像里近乎透明；
%   重建振幅图里相位粒子只剩一圈"幽灵环"（干涉伪影，Python 实跑对比度 0.26，且非真实形貌）。
% - 相位信息藏在同轴干涉条纹里：sqrt 强度 + 反向传播读出的相位被背景干涉压低
%   （实跑 0.589 rad vs 真值 1.5 rad，非线性失真）。
% - SHM（4.1）记录的电子全息图是复值的（ic/is 双通道可合成复全息图，ch3 式 3.5-5），
%   相位通道一就读出定量相位——这正是 ch5"相干模式/定量相位衬度"的入口。
%
% 差异提示：
% 1) 全书最大的坑：书的前向变换是 +j 指数，严格实现用 ifft2*N²；逆变换 fft2/N²（文件末尾封装）。
%    你旧版 Amplitude_and_Phase_Object.m 的 ifft2*N² + fft2/N² 组合正是书协议（与 2.3 节一致）。
% 2) 连通域定位用 regionprops（Image Processing Toolbox）；去 toolbox 的简化版见 Python 版。
% 3) 旧版振幅粒子取值不一致（.m 用 0.1、.py 用 0），本版统一 0.1。

clear
N = 512; lambda0 = 532e-9; dx = 5e-6; L = N*dx; z_prop = 0.05; k0 = 2*pi/lambda0;
x = ((0:N-1)-N/2)*dx;
[X, Y] = meshgrid(x, x);

% ---- 1. 物体：左 振幅粒子（0.1 吸收），右 相位粒子（1.5 rad 延迟） ----
U_obj = ones(N);
radius = 40e-6;
mask_amp = ((X + 0.4e-3).^2 + Y.^2) <= radius^2;     % 左：纯振幅（强吸收、无相位）
U_obj(mask_amp) = 0.1;
mask_ph = ((X - 0.4e-3).^2 + Y.^2) <= radius^2;      % 右：纯相位（无吸收、延迟 1.5 rad）
U_obj(mask_ph) = exp(1j*1.5);

% ---- 2. 角谱正/反向传播（书式协议，自动滤消逝波） ----
df = 1/L;
fx = ((0:N-1)-N/2)*df;
[KY, KX] = meshgrid(fx, fx);
root = 1 - (KX.^2 + KY.^2)/k0^2;  root(root < 0) = 0;   % 滤除消逝波
H_fwd = exp(-1j*k0*z_prop*sqrt(root));
H_bwd = exp(+1j*k0*z_prop*sqrt(root));                   % 反向传播：z 取 -z_prop

U_det = ift2_book(ft2_book(U_obj).*H_fwd);               % 传播到相机面
I_holo = abs(U_det).^2;                                  % 相机只记强度
U_rec = ift2_book(ft2_book(sqrt(I_holo)).*H_bwd);        % sqrt 强度当复振幅估计，反向传播

% ---- 3. 定量对比：振幅通道里"看得见"谁？相位通道里"读得出"谁？ ----
amp = abs(U_rec);
ph = angle(U_rec);
bg = mean(amp((X.^2 + Y.^2) > (1.0e-3)^2), 'all');       % 远场背景振幅
c_amp_obj = (bg - mean(amp(mask_amp), 'all'))/bg;        % 振幅粒子在重建振幅中的对比度
c_ph_obj = abs(mean(amp(mask_ph), 'all') - bg)/bg;       % 相位粒子在重建振幅中的对比度
ph_dome = angle(mean(U_rec(mask_ph), 'all'));            % 相位粒子中心读出的相位

fprintf('--- 重建通道对比（同轴 sqrt 强度重建） ---\n');
fprintf('振幅粒子（透过率 0.1）：重建振幅对比度 = %.3f（清晰黑盘，形貌真实）\n', c_amp_obj);
fprintf('相位粒子（延迟 1.5 rad，零吸收）：重建振幅对比度 = %.3f（全部来自幽灵环干涉伪影，不是形貌）\n', c_ph_obj);
fprintf('相位通道在相位粒子中心读出 = %.3f rad（真值 1.5，被背景干涉压低、非线性）\n', ph_dome);
fprintf('结论：振幅通道会误导（伪影像形貌）、相位通道读数失真 -> 定量相位成像必须走复全息图路线\n');

% ---- 4. 可视化 ----
figure('Position', [100, 100, 1300, 420]);
subplot(1,3,1); imagesc(x*1e3, x*1e3, abs(U_obj)); colormap gray; axis image; colorbar
title('1. Object amplitude (0.1 absorber | 1.5 rad dome)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(1,3,2); imagesc(x*1e3, x*1e3, amp); axis image; colorbar
title('2. Reconstructed amplitude (absorber clear, phase-ghost ring)'); xlabel('x (mm)'); ylabel('y (mm)')
subplot(1,3,3); imagesc(x*1e3, x*1e3, ph); colormap jet; axis image; colorbar
caxis([-pi pi])
title('3. Reconstructed phase (1.5 rad dome, corrupted)'); xlabel('x (mm)'); ylabel('y (mm)')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function F = ft2_book(U)
% 书式前向变换（+j 指数）
n = size(U, 1);
F = fftshift(ifft2(ifftshift(U))) * n^2;
end

function u = ift2_book(A)
% 书式逆变换（-j 指数）
n = size(A, 1);
u = fftshift(fft2(ifftshift(A))) / n^2;
end
