% ===== 书 5.3 · PSF 工程（式 5.3-1 ~ 5.3-7）：反向曲率球面波把横向分辨率提高 2 倍 =====
% 标准 OSH（5.3a）：p1=1（球面波@w0）+ p2=delta（平面波@w0+Omega），TDFZP 扫描（式 5.3-1）。
%   针孔物体的 cos 通道全息图（式 5.3-2）：i_c ~ sin[k0 r^2/(2 z0)]，重建像距 z0，
%   NA = r/z0，分辨率 lam0/NA。
% PSF 工程（5.3b）：把 p2 换成点源实现（点源 p 在光瞳面前 z_i 处 -> p2 = exp[-j k0 r^2/(2 z_i)]，
%   式 5.3-3），物面上得到会聚球面波（式 5.3-4：a = (jk0/2 pi z1) exp[+j k0 r^2/(2 z1)]）。
%   成像条件（式 5.3-5，原书印的就是这个形式，PDF 300 dpi 放大确认无误）：
%       1/(f+z_i) + 1/(f+z0+z1) = 1/f   <=>   f^2 = z_i (z0+z1)
%   选 z_i=2f、z0=z1=f/4 精确满足。波动光学复核：焦平面会聚波半径 R_f = f^2/z_i，
%   再传 z0 后半径 R_f - z0 = z1 —— 与 (5.3-5) 等价。物面上是"曲率相反的球面波对"，
%   扫描强度（式 5.3-6）= A' + B' sin[k0 r^2/z0 + Omega t] —— 啁啾率翻倍，
%   新全息图（式 5.3-7）i_c ~ sin[k0 r^2/z0]，重建像距 z0/2，NA = r/(z0/2)，分辨率 lam0/(2 NA)。
%
% 本脚本验证：a) 式 5.3-5 成像条件（书本数字 + 等价形式 f^2=z_i(z0+z1)）；b) 式 5.3-3 -> 5.3-4
%   会聚波曲率 R = f^2/z_i - z0；c) 两张全息图的聚焦距离；d) 焦斑 FWHM 之比 ~ 2（Airy 理论对照）。
%
% 严格书约定：前向变换 +j 指数 -> ifft2(...)*N^2；逆变换 -j 指数 -> fft2(...)/N^2。
% MATLAB / Python 差异提示：
% 1) 参数设计（三个约束）：角谱传递函数无混叠 d < N*dx^2/lam0；菲涅尔数 N_F=r^2/(lam0*z0)~4
%    （太小时低菲涅尔数"焦移"明显，Li & Wolf 效应）；焦斑 FWHM>=数像素可采样。
% 2) MATLAB 的 besselj(1,x) 对应 scipy.special.j1（第一类一阶贝塞尔）。

clear
N = 1024; dx = 2.0; lam0 = 0.5; k0 = 2*pi/lam0;
z0 = 7040.0; r_ap = 120.0;            % N_F = r^2/(lam0*z0) = 4.1

x = ((0:N-1) - N/2)*dx;
[X, Y] = meshgrid(x, x);
R2 = X.^2 + Y.^2;
KX = 2*pi*(0:N-1)/(N*dx); KX(KX >= pi/dx) = KX(KX >= pi/dx) - 2*pi/dx;
[KX2, KY2] = meshgrid(KX, KX);
K2 = KX2.^2 + KY2.^2;
row = N/2 + 1;                        % 中心行（MATLAB 1-based；Python 用 0-based 的 N//2）

fwd = @(u) ifft2(u) * N^2;            % 书前向变换（+j 指数）
invf = @(U) fft2(U) / N^2;            % 书逆变换（-j 指数）
propagate = @(u, d) invf( fwd(u) .* exp(-1j*k0*d) .* exp(1j*d*K2/(2*k0)) );  % 式 (2.3-13)

ap = R2 <= r_ap^2;                    % 全息图孔径（半径 r）

% ----- 1. 式 5.3-5 成像条件（原书即印 1/(f+z_i)+1/(f+z0+z1)=1/f，无误）-----
% 交叉验证：a) 书本数字代入；b) 等价形式 f^2 = z_i(z0+z1)（交叉相乘可证）；
% c) 波动光学链 R_f = f^2/z_i -> R = R_f - z0 = z1（见第 2 步），与 b) 恒等。
% （若把 p2 下标误读成 z1，则条件变成 f^2 = z1(z0+z1)，书数 z_i=2f, z0=z1=f/4 不满足——
%   故 p2 下标必为 z_i：由"点源放在光瞳面前 z_i 处"实现，见书 p.146。）
f = 8000.0; zi = 2*f; z0c = f/4; z1c = f/4;   % 书中选值：z_i=2f, z0=z1=f/4
lhs = 1/(f + zi) + 1/(f + z0c + z1c);         % 式 (5.3-5) 左端
fprintf('--- 式 5.3-5 成像条件核验（f=%d, z_i=2f, z0=z1=f/4）---\n', f);
fprintf('1/(f+z_i)+1/(f+z0+z1) = %.6f /f (应 = 1.000000，书本数字精确成立)\n', lhs*f);
fprintf('等价形式：z_i = f^2/(z0+z1) = %.4f f (应 = 2f，一致)\n', f^2/(z0c+z1c)/f);

% ----- 2. 核验式 5.3-3 -> 5.3-4：物面会聚波的曲率半径 = f^2/z_i - z0（= z1）-----
% 点源（p2 的实现）经透镜成像：焦平面场是会聚 chirp，半径 f^2/z_i；
% 再卷积 h(x,y;z0) 传到物面（式 5.3-3），应得会聚波半径 z1。
R_f = f^2/zi;                                          % 焦平面会聚波半径（本参数下 = 4000）
u_f = exp(1j*k0*R2/(2*R_f));                           % 焦平面场（会聚）
a_field = propagate(u_f, z0c);                         % 式 (5.3-3)：再卷积 h(x,y;z0)

seg = (row-60):(row+60);                               % 只在采样充分的中心区拟合
ph = unwrap(angle(a_field(row, seg)));
rho = x(seg);
cf = polyfit(rho.^2, ph, 1);                           % 相位 ~ cf(1)*r^2 + cf(2)
R_num = k0/(2*cf(1));
fprintf('--- 式 5.3-4 会聚波曲率核验（点源实现 z_i=2f）---\n');
fprintf('理论 R = f^2/z_i - z0 = %.1f  (即 z1 = %.1f)\n', f^2/zi - z0c, z1c);
fprintf('数值拟合 R_num        = %.1f  (相对误差 %.2f%%)\n', R_num, ...
    abs(R_num - (f^2/zi - z0c))/(f^2/zi - z0c)*100);

% ----- 3. 两张针孔全息图（式 5.3-2 / 5.3-7 的复全息图版，取 +chirp 支）-----
H_std = ap .* exp(1j*k0*R2/(2*z0));    % 标准 TDFZP 记录，重建像距 z0
H_eng = ap .* exp(1j*k0*R2/z0);        % 反向曲率对记录，重建像距 z0/2
s_std = ap .* sin(k0*R2/(2*z0));       % 式 (5.3-2) 的 cos 通道（显示用）
s_eng = ap .* sin(k0*R2/z0);           % 式 (5.3-7)

d_grid = linspace(0.4*z0, 1.1*z0, 29); % 聚焦距离扫描（覆盖 z0/2 与 z0，且 < AS 采样上限）
sharp_std = arrayfun(@(d) max(abs(propagate(H_std, d)), [], 'all')^2, d_grid);
sharp_eng = arrayfun(@(d) max(abs(propagate(H_eng, d)), [], 'all')^2, d_grid);
[~, i_std] = max(sharp_std);  [~, i_eng] = max(sharp_eng);
fprintf('--- 聚焦距离扫描 ---\n');
fprintf('标准全息图最锐距离 = %.0f (几何焦点 z0 = %.0f；N_F=%.1f 的低菲涅尔数焦移使峰值略前移)\n', ...
    d_grid(i_std), z0, r_ap^2/(lam0*z0));
fprintf('工程全息图最锐距离 = %.0f (几何焦点 z0/2 = %.0f)\n', d_grid(i_eng), z0/2);

% ----- 4. 焦斑 FWHM（线性插值找半高交点）-----
U_std = propagate(H_std, z0);          % 标准焦斑（z0 处）
U_eng = propagate(H_eng, z0/2);        % 工程焦斑（z0/2 处）
[w_std, w_eng] = deal_fwhm2(abs(U_std(row,:)).^2, abs(U_eng(row,:)).^2, dx);
fprintf('--- 焦斑 FWHM（横向分辨率）---\n');
fprintf('标准：FWHM = %.2f（Airy 理论 %.2f）—— 分辨率 lam0/NA，NA = r/z0\n', ...
    w_std, 3.232*z0/(k0*r_ap));
fprintf('工程：FWHM = %.2f（Airy 理论 %.2f）—— 分辨率 lam0/(2NA)\n', ...
    w_eng, 3.232*(z0/2)/(k0*r_ap));
fprintf('FWHM 之比 = %.3f（理论 2.0：超越孔径 Rayleigh 极限 2 倍）\n', w_std/w_eng);

% ----- 5. 可视化 -----
figure('Position', [40 40 1250 950]);
subplot(2, 2, 1); imagesc(x, x, s_std); axis xy equal tight; xlim([-250 250]); ylim([-250 250]);
colormap(gca, gray); title('1. Std hologram i_c ~ sin(k0 r^2 / 2 z0)  (Eq.5.3-2)');
xlabel('x'); ylabel('y');
subplot(2, 2, 2); imagesc(x, x, s_eng); axis xy equal tight; xlim([-250 250]); ylim([-250 250]);
colormap(gca, gray); title('2. Engineered i_c ~ sin(k0 r^2 / z0)  (Eq.5.3-7, chirp x2)');
xlabel('x'); ylabel('y');
subplot(2, 2, 3);
Ist = abs(U_std(row,:)).^2; Ien = abs(U_eng(row,:)).^2;
plot(x, Ist/max(Ist), x, Ien/max(Ien), 'LineWidth', 1.2); xlim([-30 30]); grid on;
legend(sprintf('std @z0, FWHM=%.1f', w_std), sprintf('eng @z0/2, FWHM=%.1f', w_eng), 'Location', 'northeast');
title('3. Focal spot cross-sections (FWHM ratio ~ 2)'); xlabel('x'); ylabel('normalized intensity');
subplot(2, 2, 4);
plot(d_grid/z0, sharp_std/max(sharp_std), 'o-', d_grid/z0, sharp_eng/max(sharp_eng), 's-', 'MarkerSize', 4);
grid on; legend('std (peak @ z0)', 'eng (peak @ z0/2)', 'Location', 'northwest');
title('4. Focus sharpness vs reconstruction distance');
xlabel('reconstruction distance d / z0'); ylabel('peak |U|^2 (norm.)');

% ===== 局部函数 =====
function w = deal_fwhm(prof, step)     % 单条剖面的半高全宽（线性插值）
    I = prof / max(prof);
    half = 0.5;
    pk = find(I == max(I), 1);
    hi = pk;
    while hi < length(I) && I(hi+1) >= half, hi = hi + 1; end
    xr = hi + (I(hi) - half)/(I(hi) - I(hi+1));
    lo = pk;
    while lo > 1 && I(lo-1) >= half, lo = lo - 1; end
    xl = lo - (I(lo) - half)/(I(lo) - I(lo-1));
    w = (xr - xl) * step;
end

function [w1, w2] = deal_fwhm2(p1, p2, step)
    w1 = deal_fwhm(p1, step);
    w2 = deal_fwhm(p2, step);
end
