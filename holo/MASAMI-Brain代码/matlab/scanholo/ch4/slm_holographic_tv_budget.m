% ===== 书 4.2 节 · 3-D 全息电视的三大瓶颈（式 4.2-6 ~ 4.2-10 + 数据率） =====
% A. 空间频率分辨率 -> 视角：
%    点源全息图 sin[(k0/2z0)(x²+y²)] 的局部空间频率 (4.2-6)：
%        f_local = (1/2π) d/dx[(k0/2z0)x²] = x/(λ0·z0)
%    全息图孔径半宽 x_max 处最大 (4.2-7)：f_max = x_max/(λ0·z0)
%    几何 (4.2-8)：NA = sin(θ/2) = x_max/z0  ->  (4.2-9)：NA = λ0·f_max = λ0·f0
%    即：SLM 分辨率 f0 直接决定视角 θ = 2·arcsin(λ0·f0)。
% B. 空间分辨率（像素数）：Nyquist 采样 fs = 2·f0；l×l 全息图需
%        N = (l·fs)² = (2l·f0)² = (2l·NA/λ0)²                (4.2-10)
% C. 数据传输：N×8bit×30帧/s；HPO（只保留水平视差）可把像素数砍到 512×(2l·NA/λ0)。
%
% 本脚本复核书中的全部数值，并验证书 4.2-C 的一处笔误：
%    书印 HPO 像素公式 512×(2l×2NA/λ0) —— 代入得 32.4 million，与书自己给的
%    "approximately 17 million" 矛盾；正确式应为 512×(2l·NA/λ0) ≈ 16.2 million（见笔记勘误）。
%
% 差异提示：纯数值脚本，无傅里叶约定问题；Python 版同目录 slm_holographic_tv_budget.py。

clear
lam0 = 0.6328e-3;         % mm（He-Ne）
l = 20.0;                 % 全息图边长 mm（20mm × 20mm）
theta_full = 60.0;        % 视角（度）
NA = sind(theta_full/2);

% ----- A. 视角 vs SLM 分辨率（Table 4.1） -----
f0_list = [8 100 500];                            % lp/mm
theta_list = 2*asind(lam0*f0_list);
book_theta = {'0.6', '6.8', '34'};
book_dev = {'EBSLM/Hamamatsu', 'PALSLM/Hamamatsu', 'not available'};
fprintf('--- A. 视角（Table 4.1，λ0=0.6328 µm） ---\n');
for ii = 1:3
    fprintf('    f0 = %5.0f lp/mm -> θ = %5.2f deg（书 Table 4.1：%s°，%s）\n', ...
            f0_list(ii), theta_list(ii), book_theta{ii}, book_dev{ii});
end
fprintf('    离轴全息 45° 记录角 -> 载频 ~1000 lp/mm：远超现有 SLM（同轴 OSH 才上得了 SLM）\n');

% ----- B. 像素数（4.2-10） -----
f0_req = NA/lam0;                       % 满足 60° 视角所需的 SLM 分辨率（4.2-9 反推）
fs = 2*f0_req;                          % Nyquist：fs = 2 f0
N_pix = (l*fs)^2;
fprintf('--- B. 像素数（式 4.2-10，l=20mm，θ=60°） ---\n');
fprintf('    所需 f0 = NA/λ0 = %.0f lp/mm，fs = 2f0 = %.0f lp/mm\n', f0_req, fs);
fprintf('    N = (l·fs)² = (2l·NA/λ0)² = %.2f billion（书：about 1.1 billion，按 λ0=0.6µm 计恰为 1.1）\n', N_pix/1e9);
fprintf('    对照 Canon D60：3072×2048 ≈ 6.3 MP（67.7 lp/mm，7.4 µm 像素）——差 160 倍\n');

% ----- C. 数据率与 HPO -----
N_book = 1.1e9;                                   % 书采用 1.1 billion
fprintf('--- C. 数据率 ---\n');
fprintf('    全视差：%.3e samples × 8 bit × 30 fps = %.2f Tbit/s（书用 1.1 billion -> 0.26 Tbit/s）\n', ...
        N_pix, N_pix*8*30/1e12);
n_hpo_book = 512*(2*l*2*NA/lam0);                 % 书印公式（多一个因子 2）
n_hpo_fix  = 512*(2*l*NA/lam0);                   % 与 17 million 自洽的式
fprintf('    HPO（512 条垂直线）：\n');
fprintf('      书印式 512×(2l×2NA/λ0) = %.1f million  <- 与书自己的 "17 million" 矛盾（笔误）\n', n_hpo_book/1e6);
fprintf('      修正式 512×(2l·NA/λ0)  = %.1f million ≈ 17 million ✓  -> %.1f Gbit/s（书：4 Gbit/s）\n', ...
        n_hpo_fix/1e6, n_hpo_fix*8*30/1e9);
fprintf('    光纤 40 Gbit/s 可达 -> 实时全息电视可行，前提是 HPO 电子全息记录（HPO-OSH，1-D TDFZP 扫描）\n');

% ----- 可视化（对应笔记图 sim-slm-budget.png） -----
figure(1)
subplot(1,2,1)
f0s = linspace(1, 600, 400);
th = 2*sind(lam0*f0s);
plot(f0s, th, 'b-', 'LineWidth', 1.5); hold on; grid on
plot(f0_list, theta_list, 'ro', 'MarkerSize', 6, 'LineWidth', 1.2)
yline(60, 'r--', '60° target');
xlabel('SLM resolution f_0 (lp/mm)'); ylabel('viewing angle θ (deg)')
title('Viewing angle vs SLM resolution (Eq. 4.2-9)')

subplot(1,2,2)
rates = [N_book*8*30/1e9, n_hpo_fix*8*30/1e9];
bar(rates, 0.5); set(gca, 'YScale', 'log'); grid on
set(gca, 'XTickLabel', {'full parallax', 'HPO (512 lines)'})
ylabel('serial data rate (Gbit/s)'); ylim([1 1e4])
text(1, rates(1)*1.2, sprintf('%.0f Gbit/s', rates(1)), 'HorizontalAlignment', 'center')
text(2, rates(2)*1.2, sprintf('%.0f Gbit/s', rates(2)), 'HorizontalAlignment', 'center')
yline(40, 'k:', 'fiber 40 Gbit/s');
title('Data rate: 8 bit × 30 fps')
