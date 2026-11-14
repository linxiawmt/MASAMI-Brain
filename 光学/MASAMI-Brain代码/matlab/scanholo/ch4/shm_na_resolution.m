% ===== 书 4.1 节 · 分辨率-景深度量衡（式 4.1-1 ~ 4.1-5）与 Schilling 荧光实验的数值账 =====
% 量子力学不确定关系（书用 h 而非 hbar）：
%   横向：Δr·Δp_r >= h,  Δp_r = p0·sin(θ/2),  p0 = h/λ0
%       -> Δr >= λ0/sin(θ/2) = λ/(n0 sin(θ/2)) = λ/NA          (4.1-2)
%   纵向：Δz·Δp_z >= h,  Δp_z = p0·[1-cos(θ/2)]
%       -> Δz >= λ0/[1-cos(θ/2)] ≈ 2λ0/sin²(θ/2) = 2n0λ/NA²    (4.1-4)
%   两式结合：(Δr)²/Δz >= λ0/2                                  (4.1-5，标度律)
% 分辨率提高 2 倍 -> 景深必须缩短 4 倍。
%
% 实验参数（书 4.1 节，Schilling et al. Optics Lett. 22, 1506 (1997)）：
%   f = 150 mm（聚焦透镜 L1 焦距），D = 10 mm（被聚焦的平面波直径），λ0 = 514.5 nm（氩离子）
% 书印数值：NA ≈ 0.033，Δr ≈ 18.5 µm，Δz ≈ 1028.4 µm。
% 实跑重算：NA = 0.0333 -> Δr = 15.4 µm、Δz = 0.93 mm（书值对应 NA≈0.028/0.0316，舍入不一致）。
%
% 差异提示：纯数值脚本，无傅里叶约定问题；Python 版同目录 shm_na_resolution.py。

clear
f1  = 150e-3;        % L1 焦距 (m)
D   = 10e-3;         % 平面波直径 (m)
lam = 514.5e-9;      % 氩离子 514.5 nm 线（真空中）
n0  = 1.0;           % 样本空间近似空气/水

half_angle = atan((D/2)/f1);          % θ/2
NA = n0 * sin(half_angle);
d_r = lam / NA;                       % (4.1-2) 横向分辨率
d_z = 2 * n0 * lam / NA^2;            % (4.1-4) 景深

fprintf('θ/2 = %.4f deg,  NA = %.4f\n', rad2deg(half_angle), NA);
fprintf('式(4.1-2): Δr = λ0/NA = %.1f µm   （书印 18.5 µm）\n', d_r*1e6);
fprintf('式(4.1-4): Δz = 2λ0/NA² = %.1f µm = %.2f mm（书印 1028.4 µm）\n', d_z*1e6, d_z*1e3);
fprintf('15 µm 荧光小球 vs Δr = %.1f µm：接近横向极限，可分辨\n', d_r*1e6);
fprintf('式(4.1-5) 标度律：(Δr)²/Δz = %.4f µm ≈ λ0/2 = %.4f µm\n', d_r^2/d_z*1e6, lam/2*1e6);

% ---- (4.1-5) 标度律曲线：Δz(Δr) 双曲线 ----
r = linspace(5e-6, 60e-6, 400);
z_hyp = 2 * r.^2 / lam;               % 由 (Δr)²/Δz = λ0/2 反解

figure(1)
subplot(1,2,1)
plot(r*1e6, z_hyp*1e3, 'b-', 'LineWidth', 1.5); hold on
plot(d_r*1e6, d_z*1e3, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)
grid on
xlabel('lateral resolution Δr (µm)'); ylabel('depth of focus Δz (mm)')
title('Uncertainty trade-off in microscopy (Eq. 4.1-5)')
legend('(\Deltar)^2/\Deltaz = \lambda_0/2 (4.1-5)', ...
       sprintf('Schilling setup (%.1f µm, %.2f mm)', d_r*1e6, d_z*1e3), 'Location', 'northwest')

% ---- TDFZP 条纹间距账：扫描场内最细条纹（决定采样要求） ----
subplot(1,2,2)
xv = linspace(0, D/2, 300);                            % 扫描场半径（受平面波直径限制）
fringe_period = lam * f1 ./ xv;                        % 局部周期 λ0·f1/x
semilogy(xv*1e3, fringe_period*1e6, 'g-', 'LineWidth', 1.5); hold on
yline(2*d_r*1e6, 'r--', '2·Δr', 'LineWidth', 1.2);
grid on
xlabel('radius in scan field (mm)'); ylabel('local FZP fringe period (µm)')
title('TDFZP fringe period across scan field')
