%% ===== 书 Schmidt PM199 第 1 章 =====
%  单缝衍射 / 双缝干涉：菲涅耳近场 + 夫琅禾费远场
%  对应原书式 (1.57)（菲涅耳衍射积分）与式 (1.67)、(1.71)（夫琅禾费衍射）。
%  重写自用户旧稿 InterferenceDiffraction.m / 双缝干涉和单缝衍射.py：
%    - 旧稿 propFresnelTF 里漏写了 fft2（H 根本没乘到频谱上），本版已修正；
%    - 统一采用"胞元中心"采样：x 偏移半像素，40 um 缝恰占 20 个胞元，
%      离散谱零点严格落在连续理论位置 m*lambda*z/Dx 上。
%  MATLAB 注意：数组从 1 起；fftshift 把 DC 移到第 N/2+1 个元素。
%  需要 R2016b+（脚本内含本地函数）。
clear; clc; close all;

%% 1. 物理与网格参数
N    = 1024;            % 网格点数（N x N）
dx   = 2e-6;            % 采样间隔 2 um
wvl  = 532e-9;          % 波长 532 nm
k    = 2*pi/wvl;

x1 = ((1:N) - N/2 - 0.5) * dx;            % 胞元中心坐标（与 Python (arange(N)-N//2+0.5)*dx 相同）
[X1, Y1] = meshgrid(x1, x1);              % MATLAB 与 numpy 的 meshgrid 同名同义

% 单缝：宽 Dx、高 Dy；双缝：中心距 dx_slit（> Dx，两缝不重叠）
Dx      = 40e-6;
Dy      = 400e-6;
dx_slit = 120e-6;

U_single = double(abs(X1) <= Dx/2 & abs(Y1) <= Dy/2);      % rect(x/Dx)rect(y/Dy)
U_double = double(abs(abs(X1) - dx_slit/2) <= Dx/2 & abs(Y1) <= Dy/2);

z_near = 5e-3;          % 近场观察距离 5 mm（菲涅耳区）
z_far  = 80e-3;         % 远场观察距离 80 mm（夫琅禾费区）

%% 2. 计算
% --- A. 菲涅耳近场（传递函数法 TF）---
I_single_fres = abs(propFresnelTF(U_single, dx, wvl, k, z_near)).^2;
I_double_fres = abs(propFresnelTF(U_double, dx, wvl, k, z_near)).^2;

% --- B. 夫琅禾费远场（式 1.67：FT + 坐标缩放 x2 = lambda*z*fx）---
[U_s_far, X2, Y2] = fraunhoferFFT(U_single, dx, wvl, z_far);
U_d_far = fraunhoferFFT(U_double, dx, wvl, z_far);
I_single_far = abs(U_s_far).^2;
I_double_far = abs(U_d_far).^2;

x1_um = x1*1e6;
x2_um = X2(1,:)*1e6;

%% 3. 数值验证（4 倍补零，y2 = 0 剖面）
[I_pad_d, x2_pad] = paddedFarProfile(U_double, dx, wvl, z_far);
[I_pad_s, ~     ] = paddedFarProfile(U_single, dx, wvl, z_far);

fprintf('===== 夫琅禾费衍射数值验证（4 倍补零，y2 = 0 剖面）=====\n');

% (1) 左右对称性：在镜像点 -x2 处插值取值再比较。
% 注意：偶数 N 的中心化频率轴是 -N/2 ... N/2-1 个 bin，左右差一个 bin，
% 直接 I 与 flip(I) 比较会错位一个 bin，必须插值到 -x2。
prof_d = I_double_far(N/2+1, :);                       % y=0 行（MATLAB 索引 N/2+1）
mirror = interp1(x2_um, prof_d, -x2_um, 'linear', 'extrap');
sym = max(abs(prof_d - mirror))/max(prof_d);
fprintf('(1) 双缝远场左右对称性 max|I(x2)-I(-x2)|/max(I) = %.2e\n', sym);

% (2) 单缝暗纹：sinc 零点（严格）理论位置 x2 = m*lambda*z/Dx
dark_th = (1:7) * wvl*z_far/Dx * 1e6;
win = x2_pad > 0 & x2_pad < dark_th(end)*1.1;
xm = x2_pad(win); pm = I_pad_s(win);
idx = find(pm(2:end-1) < pm(1:end-2) & pm(2:end-1) < pm(3:end)) + 1;
nb2 = min(numel(dark_th), numel(idx));
dark_num = xm(idx(1:nb2));
rel_dark = max(abs(dark_num - dark_th(1:nb2))./dark_th(1:nb2));
fprintf('(2) 单缝暗纹理论 (um): %s\n    数值极小   (um): %s  最大相对偏差 %.2e\n', ...
    mat2str(round(dark_th,1)), mat2str(round(dark_num,1)), rel_dark);

% (3) 双缝：暗纹零点严格位于 (n+1/2)*lambda*z/dx_slit 与包络零点 m*lambda*z/Dx；
%     亮纹极大受单缝包络斜坡影响，会从 cos^2 极大 n*lambda*z/dx_slit 向中心偏移
%     （第一亮纹约 341 um 而非 354.7 um），故与式 (1.71) 解析剖面的极大对比。
sp = wvl*z_far/dx_slit*1e6;                            % cos^2 条纹间距 (um)
cos_zero_th = ((0:5)+0.5) * sp;
env_zero_th = (1:6) * wvl*z_far/Dx * 1e6;
sub = x2_pad > 0 & x2_pad < 6.5*sp;
xd = x2_pad(sub); pd = I_pad_d(sub);
frac = zeros(5,1);
zlist = [cos_zero_th(1:3), env_zero_th(1:2)];
for jj = 1:5
    [~, i0] = min(abs(xd - zlist(jj)));
    inw = abs(xd - zlist(jj)) < 0.9*sp;
    frac(jj) = pd(i0)/max(pd(inw));
end
fprintf('(3a) 双缝暗纹深度：理论零点处 I/邻近条纹峰 = %s（全部 << 1）\n', ...
    mat2str(round(frac, 5)));

x2_fine = 0 : 0.2 : 6.2*sp;                            % 解析式细网格 (um)
pa_fine = cos(pi*dx_slit*(x2_fine*1e-6)/(wvl*z_far)).^2 ...
        .* mysinc(Dx*(x2_fine*1e-6)/(wvl*z_far)).^2;   % 式(1.71) 沿 y2=0（去掉常数）
idx_b = find(pa_fine(2:end-1) > pa_fine(1:end-2) & pa_fine(2:end-1) > pa_fine(3:end)) + 1;
bright_ana = x2_fine(idx_b);
idx_n = find(pd(2:end-1) > pd(1:end-2) & pd(2:end-1) > pd(3:end)) + 1;
bright_num = xd(idx_n);
nb = min(numel(bright_num), numel(bright_ana));
rel_bright = max(abs(bright_num(1:nb) - bright_ana(1:nb))./bright_ana(1:nb));
fprintf('(3b) 亮纹 cos^2 极大（朴素）(um): %s\n     解析式(1.71)极大 (um): %s  <- 被包络拉向中心\n     数值极大 (um): %s  最大相对偏差 %.2e\n', ...
    mat2str(round((1:6)*sp,1)), mat2str(round(bright_ana(1:nb),1)), ...
    mat2str(round(bright_num(1:nb),1)), rel_bright);
cand = [3 6];
miss = cand(arrayfun(@(n) min(abs(bright_num - n*sp)) > 30, cand));
fprintf('     缺级检验：缺失级次 = %s（理论预期 dx_slit/Dx=3 -> [3 6]）\n', mat2str(miss));

% (4) 与式 (1.71) 解析解对比：|U| 沿 y2 = 0，取显示窗口 ±800 um（未补零原网格）
U_ana = 2*Dx*Dy/(wvl*z_far) * cos(pi*dx_slit*X2/(wvl*z_far)) ...
      .* mysinc(Dx*X2/(wvl*z_far)) .* mysinc(Dy*Y2/(wvl*z_far));
win4 = abs(x2_um) <= 800;
err = max(abs(abs(U_d_far(N/2+1, win4)) - abs(U_ana(N/2+1, win4)))) / max(abs(U_d_far(N/2+1,:)));
fprintf('(4) FFT vs 式(1.71)（±800 um 内）max||U|-|U_ana||/max|U| = %.2e\n', err);

%% 4. 绘图（caxis 截到 2%：压掉 0 级亮斑，突出条纹）
figure('Name','单缝/双缝衍射','Position',[60 60 1150 850]);
panel = {I_single_fres, x1_um, '单缝·菲涅耳近场 z=5 mm'; ...
         I_single_far,  x2_um, '单缝·夫琅禾费远场 z=80 mm'; ...
         I_double_fres, x1_um, '双缝·菲涅耳近场 z=5 mm'; ...
         I_double_far,  x2_um, '双缝·夫琅禾费干涉 z=80 mm'};
for jj = 1:4
    subplot(2,2,jj);
    imagesc(panel{jj,2}, panel{jj,2}, panel{jj,1}); axis xy; axis image;
    caxis([0, 0.02*max(panel{jj,1}(:))]); colormap gray; colorbar;
    xlim([-800 800]); ylim([-800 800]);
    title(panel{jj,3}); xlabel('x (\mum)'); ylabel('y (\mum)');
end

figure('Name','双缝夫琅禾费剖面','Position',[80 80 900 420]);
plot(x2_pad, I_pad_d/max(I_pad_d), 'b-', 'LineWidth', 1.2); hold on;
plot(x2_pad, mysinc(Dx*(x2_pad*1e-6)/(wvl*z_far)).^2, 'g--', 'LineWidth', 1.2);
xlim([-800 800]); ylim([-0.02 1.05]); grid on;
xlabel('远场坐标 x_2 (\mum)'); ylabel('归一化光强');
legend('双缝 FFT 数值解（4 倍补零）', '单缝包络 sinc^2(D_x x_2/\lambdaz)', 'Location', 'northeast');
title('双缝夫琅禾费干涉：\Deltax = 3D_x → ±3 级缺级');

%% ===== 本地函数 =====
function U_out = propFresnelTF(U_in, delta, wavelength, k, z)
    % 菲涅耳衍射·传递函数法：U2 = IFT2{ FT2{U1} .* H }，
    % H = exp(ikz) exp(-i*pi*lambda*z*(fx^2+fy^2))。
    % 适用条件：z < 4*dx^2/lambda（本例 30 mm）。
    [M, N] = size(U_in);
    fx = ((1:N) - N/2 - 1) / (N*delta);      % 与 Python (arange(N)-N//2)/(N*delta) 相同
    fy = ((1:M) - M/2 - 1) / (M*delta);
    [FX, FY] = meshgrid(fx, fy);
    H = exp(1i*k*z) * exp(-1i*pi*wavelength*z*(FX.^2 + FY.^2));
    U_ft = fftshift(fft2(fftshift(U_in)));   % 中心化谱
    U_out = ifftshift(ifft2(ifftshift(U_ft .* H)));
end

function [U_out, X2, Y2] = fraunhoferFFT(U_in, delta, wavelength, z)
    % 夫琅禾费衍射（式 1.67）：观察面坐标 x2 = lambda*z*fx。
    [M, N] = size(U_in);
    fx = ((1:N) - N/2 - 1) / (N*delta);
    fy = ((1:M) - M/2 - 1) / (M*delta);
    [FX, FY] = meshgrid(fx, fy);
    x2 = wavelength*z*FX;  y2 = wavelength*z*FY;
    G = fftshift(fft2(fftshift(U_in))) * delta^2;      % 连续 FT 的离散近似
    U_out = exp(1i*2*pi/wavelength*z) * exp(1i*2*pi/wavelength*(x2.^2+y2.^2)/(2*z)) ...
          / (1i*wavelength*z) * G;
    X2 = x2; Y2 = y2;
end

function [I_pad, x2_pad] = paddedFarProfile(U_in, delta, wavelength, z)
    % 4 倍补零后再取远场 |U| 剖面：极值定位精度提高约 4 倍。
    % 补零在 fftshift 后的一维序列上进行：正坐标样本（前半段）放头部，
    % 负坐标样本（后半段）接到尾部，保持采样列等距连续、相位参考点不动。
    N = length(U_in); Np = 4*N;
    h = fftshift(U_in(N/2+1, :));            % y=0 行
    hp = zeros(1, Np);
    hp(1:N/2) = h(1:N/2);
    hp(Np-N/2+1:end) = h(N/2+1:end);
    F = fft(hp) * delta;
    fk = (0:Np-1)/(Np*delta);
    fk(fk > 1/(2*delta)) = fk(fk > 1/(2*delta)) - 1/delta;   % 折回负频率
    I_pad = abs(F).^2;
    x2_pad = fk * wavelength * z * 1e6;
end

function s = mysinc(x)
    % 归一化 sinc：sinc(x) = sin(pi x)/(pi x)，与原书 Appendix A 定义一致。
    % （不直接用 MATLAB 内置 sinc，避免依赖 Signal Processing Toolbox）
    s = ones(size(x));
    nz = x ~= 0;
    s(nz) = sin(pi*x(nz))./(pi*x(nz));
end
