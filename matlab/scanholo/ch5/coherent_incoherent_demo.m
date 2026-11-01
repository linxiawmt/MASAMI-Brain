% ===== 书 5.1 · 相干模式 vs 非相干模式记录（式 5.1-6b / 5.1-6c，例 5.1 式 5.1-9）=====
% 演示广义双瞳处理器的核心结论（Poon & Indebetouw 2003）：
%   - 非相干模式（开掩膜 M=1）：记录的是强度物体 |T|^2，纯相位结构完全不可见（式 5.1-6c）；
%   - 相干模式（针孔掩膜 M=delta 且 p1=delta）：记录的是振幅物体 T 本身（式 5.1-6b），
%     p2=1 时就是例 5.1 的"复菲涅耳波带片全息图"（式 5.1-9）。
% 物体设计：三缝（振幅物体）+ 纯相位圆盘（|T|=1，相位 1.5 rad）。
%
% 严格书约定：前向变换 +j 指数 -> MATLAB 用 ifft2(...)*N^2 实现（ifft2 自带 1/N^2，
%   乘回 N^2 恰好是无归一化的 +j 求和）；逆变换 -j 指数 -> fft2(...)/N^2。
%   （numpy 的 fft2/ifft2 归一化与 MATLAB 相同，故 Python 版用同样的 ifft2*N^2 / fft2/N^2。）
%
% MATLAB / Python 差异提示：
% 1) meshgrid 与 numpy.meshgrid 默认 'xy' 一致：数组都是 (行=y, 列=x)。
% 2) imshow 坐标原点：MATLAB imagesc 默认 y 向上；Python plt.imshow 需 origin='lower'。
% 3) 相位显示里振幅为零处置 NaN（显示为白），避免 angle(0±eps) 的随机噪声误导。

clear
N = 512; dx = 1.0; lam0 = 0.5; k0 = 2*pi/lam0; z0 = 2000.0;   % 任意单位（无量纲演示）

x = ((0:N-1) - N/2)*dx;
[X, Y] = meshgrid(x, x);                                      % (行=y, 列=x)
KX = 2*pi*(0:N-1)/(N*dx);                                     % 频域网格（同 numpy.fft.fftfreq 的排列）
KX(KX >= pi/dx) = KX(KX >= pi/dx) - 2*pi/dx;                  % 后半段折回负频率
[KX2, KY2] = meshgrid(KX, KX);
K2 = KX2.^2 + KY2.^2;

% ----- 1. 物体 T(x,y;z=0)：振幅三缝 + 纯相位圆盘 -----
T = zeros(N, N);
for xc = [-180, -140, -100]                                   % 三条振幅缝（宽 12）
    T(abs(X - xc) <= 6 & abs(Y) <= 120) = 1.0;
end
disc = (X - 70).^2 + Y.^2 <= 60^2;                            % 纯相位圆盘：零吸收
T(disc) = exp(1j*1.5);                                        % 相位延迟 1.5 rad

% ----- 2. 相干模式记录（式 5.1-6b + 例 5.1 的 p2=1 -> 式 5.1-9）-----
% i_Omegap(x,y) = exp[-j k0 (z+z0)] * j k0 / (2 pi (z+z0)) * [ chirp(z+z0) * T ]
% chirp(Z) = exp[-j k0 (x^2+y^2)/(2Z)]，其频谱是 exp[+j Z k^2 / (2 k0)]（书前向 +j）。
A_coh = 1j*k0/(2*pi*z0) * exp(-1j*k0*z0);                     % 式 (5.1-9) 前因子
Coh = A_coh .* ifft2( fft2(T) .* exp(1j*z0*K2/(2*k0)) ) * N^2;  % 复 FZP 全息图

% 相干重建：全息图反向共轭传播，乘回前因子
Rc = fft2( ifft2(Coh) .* exp(-1j*z0*K2/(2*k0)) .* N^2 ) / N^2 / A_coh;

% ----- 3. 非相干模式记录（式 5.1-6c，M=1；p1=1 球面波 + p2=delta 平面波，同 ch3 标准 OSH）-----
% 取 ch3 (3.5-5) 的 +chirp 分支：H = |T|^2 * exp[+j k0 r^2/(2 z0)] * k0/(2 pi z0)，
% 该核的频谱是 exp[-j z0 k^2 / (2 k0)]。
I2 = abs(T).^2;
A_inc = k0/(2*pi*z0);
Inc = A_inc .* ifft2( fft2(I2) .* exp(-1j*z0*K2/(2*k0)) ) * N^2;  % 复 FZP 全息图（强度物体）

Ri = fft2( ifft2(Inc) .* exp(1j*z0*K2/(2*k0)) .* N^2 ) / N^2 / A_inc;  % 重建回 |T|^2

% ----- 4. 数值验证 -----
mycorr = @(a, b) sum((a(:)-mean(a(:))).*(b(:)-mean(b(:)))) / ...
    norm(a(:)-mean(a(:))) / norm(b(:)-mean(b(:)));            % 皮尔逊相关（免工具箱）
c_coh  = mycorr(abs(Rc), abs(T));                             % 相干重建振幅 vs |T|
c_inc  = mycorr(abs(Ri), I2);                                 % 非相干重建 vs |T|^2
ph_coh = angle(mean(Rc(disc)));                               % 圆盘上读出的相位
ph_inc = angle(mean(Ri(disc)));
fprintf('--- 式 5.1-6b（相干） vs 5.1-6c（非相干）记录对比 ---\n');
fprintf('相干重建：corr(|Rc|, |T|)    = %+.4f\n', c_coh);
fprintf('相干重建：圆盘相位读出      = %+.4f rad（真值 1.5）\n', ph_coh);
fprintf('非相干重建：corr(|Ri|,|T|^2) = %+.4f\n', c_inc);
fprintf('非相干重建：圆盘相位读出    = %+.4f rad（只记录了 |T|^2，相位丢失）\n', ph_inc);
disp('结论：掩膜从"开"换成"针孔"，记录对象从 |T|^2 变成 T —— 3-D 复相干处理的入口');

% ----- 5. 可视化 -----
phT = angle(T);  phT(abs(T) < 1e-9) = NaN;                    % 振幅为零处不显示相位
phRc = angle(Rc); phRc(abs(Rc) < 1e-9) = NaN;
phRi = angle(Ri); phRi(abs(Ri) < 1e-9) = NaN;

figure('Position', [50 50 1500 860]);
subplot(2, 3, 1); imagesc(x, x, phT); axis xy equal tight; colormap(gca, jet); colorbar;
title('1. Object arg(T): slits (amp) + disc (1.5 rad, |T|=1)'); xlabel('x'); ylabel('y');
subplot(2, 3, 2); imagesc(x, x, abs(Coh)); axis xy equal tight; colormap(gca, gray);
title('2. Coherent hologram |i_Omegap| (Eq.5.1-9, complex FZP)'); xlabel('x'); ylabel('y');
subplot(2, 3, 3); imagesc(x, x, phRc); axis xy equal tight; colormap(gca, jet); colorbar;
title('3. Coherent recon: arg(Rc) reads 1.5 rad on disc'); xlabel('x'); ylabel('y');
subplot(2, 3, 4); imagesc(x, x, I2); axis xy equal tight; colormap(gca, gray); caxis([0 1.2]);
title('4. Intensity object |T|^2 (what M=1 mode sees)'); xlabel('x'); ylabel('y');
subplot(2, 3, 5); imagesc(x, x, abs(Ri)); axis xy equal tight; colormap(gca, gray);
title('5. Incoherent recon: |Ri| ~ |T|^2 (slits + flat disc)'); xlabel('x'); ylabel('y');
subplot(2, 3, 6); imagesc(x, x, phRi); axis xy equal tight; colormap(gca, jet); colorbar;
title('6. Incoherent recon: arg(Ri) ~ 0 (no phase info)'); xlabel('x'); ylabel('y');
