%% ===== 书 Schmidt PM199 第 9 章（工程 1）=====
%  通过大气湍流的传播：分步束传播法 + 相位屏 + 采样 + 验证
%  （原书 9.1-9.6 节，Listings 9.1-9.8）
%  与 python/waveprop/ch9/ch9_turbulence.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准（rng 播种相同）。
%
%  原书要点：
%    分步束传播法（式 9.1-9.3）：真空衍射算符 R + 折射算符
%      T[zi,zi+1] = exp(+i*phi)（式 9.2 印成 exp(-i*phi)，按薄屏物理与
%      Listing 9.8 的 exp(+i*phz) 实现应为 + 号 —— 零均值随机屏无统计影响），
%      phi(r) = k*∫δn dz。大气 n≈1，只需在式 (8.18) 的每个部分传播平面乘屏。
%    Kolmogorov 统计（9.2.1）：Dv ~ r^(2/3)（式 9.5）、Dn（式 9.13）、
%      折射率 PSD（式 9.15-9.18）、冻结湍流（式 9.19）。
%    Rytov 微扰（9.2.2-3）：U = U0*exp(psi)（式 9.25-9.27）；平面波
%      Kolmogorov 结构函数 D_phi = 6.88 (r/r0)^(5/3)（式 9.44），
%      r0,pw = [0.423 k^2 Cn2 Dz]^(-3/5)（式 9.42）。
%    分层模型（9.2.4-5）：n 块薄屏矩匹配（式 9.65），r0i（式 9.70），
%      矩阵方程（式 9.75）+ 约束优化（fmincon）解各屏 r0i。
%    蒙特卡洛相位屏（9.3）：FT 法（Listing 9.2）+ Lane 次谐波（Listing 9.3）。
%    采样约束（9.4）：Coy 湍流展宽 D1' = D1 + c*lambda*Dz/r0,rev、
%      D2' = D2 + c*lambda*Dz/r0（式 9.84/9.85，c=2 约 97% 光），
%      约束 1/2 用 D1'/D2'，约束 3 不变；Dz_max = min(d1,dn)^2*N/lambda
%      （式 9.89）、nmin = ceil(Dz/Dz_max)+1（式 9.90）。
%
%  原书勘误（详见笔记；本脚本按修正实现）：
%    (a) 式 (9.2) T = exp(-i*phi) 应为 exp(+i*phi)（时间约定 exp(-i2*pi*nu*t)）；
%    (b) 式 (9.34) 的核应为 J1 形式（sinc 核给 0.568，J1 核给 0.49 = 式 9.49）；
%    (c) 式 (9.48) 系数应为 2*pi*k^2*Dz（书印 2*pi^2，差 pi 倍）；
%    (d) 正文"重复 40 次" vs Listing 9.8 nreals = 20（本脚本取 40）；
%    (e) Listing 9.8 第 31 行 "...." 应为 "..."。
%
%  Python 实跑结果（本脚本按相同参数/播种复现）：
%    验证 0：r0,pw = 7.03 cm、r0,sw = 12.7 cm（书正文）、sigma_chi^2,sw = 0.436；
%    验证 1：屏幕 r0i 约束优化回代 r0,sw/Rytov 相对差 < 1e-6；
%    验证 2：约束 2 给 log2 N = 8.47 -> N = 512；Dz_max = 51.2 km -> nmin = 2；
%    验证 3：真空 ROI 辐照度 std/mean ~ 1e-2 量级、准直后相位 RMS ~ 1e-2 rad；
%    验证 5/6：末屏结构函数、相干因子与理论吻合（详见笔记）。
%
%  MATLAB/Python 差异（对应 .py 内注释）：
%    - MATLAB 1-based：居中排列的零坐标样本在 N/2+1（Python N//2）；
%    - MATLAB randn 用 rng(9) 播种对应 numpy default_rng(9)；
%    - 约束优化用 fmincon（需 Optimization Toolbox）；
%    - 书 Listing 9.1 的 nsq.^8/w^16 大数幂改写为 exp(-(nsq/w^2)^8)。
clear; clc; close all;
rng(9);                                  % 对应 numpy default_rng(9)

outdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
    'docs', 'images', 'holo', 'waveprop');
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ================= 验证 0/1：几何/大气参数与屏幕优化（书 Listing 9.5） =================
fprintf('===== 验证 0：大气参数总账（书 9.5.1 节数字）=====\n');
D2 = 0.5;            % diam of the observation aperture [m]
wvl = 1e-6;          % optical wavelength [m]
k = 2*pi/wvl;        % optical wavenumber [rad/m]
Dz = 50e3;           % propagation distance [m]
DROI = 4*D2;         % diam of obs-plane region of interest [m]
D1 = wvl*Dz/DROI;    % width of central lobe of the model pt source [m]
R = Dz;              % wavefront radius of curvature [m]
Cn2 = 1e-16;         % structure parameter [m^-2/3], constant
r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);    % 式 (9.43)，常数 Cn2
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);        % 式 (9.42)
p = linspace(0, Dz, 1e3);
rytov = 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6))*(p(2)-p(1));  % 式 (9.64)
fprintf('    D1 = %.4f m, r0,pw = %.3f cm, r0,sw = %.3f cm（书 12.7 cm）\n', ...
    D1, r0pw*100, r0sw*100);
fprintf('    sigma_chi^2,sw = %.3f（书 0.436）; r0,sw/r0,pw = %.4f\n', rytov, r0sw/r0pw);

% ---- 相位屏 r0i 的约束优化（式 9.65-9.75，书 25-41 行）----
nscr = 11;
alpha = (0:nscr-1)/(nscr-1);
A = zeros(2, nscr);
A(1,:) = alpha.^(5/3);                           % r0,sw 行（式 9.72）
A(2,:) = (1-alpha).^(5/6).*alpha.^(5/6);         % Rytov 行（式 9.74）
b = [r0sw^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
x0 = (nscr/3*r0sw*ones(nscr,1)).^(-5/3);         % 初值
rmax = 0.1;                                      % 每屏 Rytov 贡献上限
x2 = rmax/1.33*(k/Dz)^(5/6)./A(2,:);             % 上界（式 9.74 反解）
x2(A(2,:)==0) = 50^(-5/3);                       % 端点屏幕名义上界 r0=50 m
fun = @(X) sum((A*X(:)-b).^2);
[Xs, fval] = fmincon(fun, x0, [], [], [], [], zeros(nscr,1), x2);
r0scrn = Xs.^(-3/5);
r0scrn(isinf(r0scrn)) = 1e6;
bp = A*Xs(:);
fprintf('===== 验证 1：屏幕 r0i 约束优化 =====\n');
fprintf('    ||A*X-b||^2 = %.3e；r0scrn [m] = %s\n', fval, mat2str(round(r0scrn,3)));
fprintf('    回代 r0,sw = %.4f cm（目标 %.4f）、Rytov = %.4f（目标 %.4f）\n', ...
    bp(1)^(-3/5)*100, r0sw*100, bp(2)*1.33*(Dz/k)^(5/6), rytov);

%% ================= 验证 2：采样分析（书 Listing 9.6） =================
fprintf('===== 验证 2：采样分析（c = 2）=====\n');
c = 2;
D1p = D1 + c*wvl*Dz/r0sw;
D2p = D2 + c*wvl*Dz/r0sw;
d1 = 10e-3;  dn = 10e-3;  N = 512;
Nneed = (wvl*Dz + D1p*dn + D2p*d1)/(2*d1*dn);    % 约束 2（式 9.87）
zmax = min([d1 dn])^2*N/wvl;                     % 式 (9.89)
nmin = ceil(Dz/zmax) + 1;                        % 式 (9.90)
fprintf('    D1p = %.4f m, D2p = %.4f m\n', D1p, D2p);
fprintf('    约束 2：N >= %.1f（log2 = %.2f，书 more than 2^8.47）\n', Nneed, log2(Nneed));
fprintf('    约束 4：Dz_max = %.1f km -> nmin = %d 平面（实际用 11 平面）\n', zmax/1e3, nmin);

%% ================= 验证 3：真空仿真（书 Listing 9.7） =================
fprintf('===== 验证 3：真空仿真 =====\n');
n = nscr;
[x1 y1] = meshgrid((-N/2:N/2-1)*d1);
[theta1 r1] = cart2pol(x1, y1);
pt = exp(-1i*k/(2*R)*r1.^2)/D1^2.*sinc(x1/D1).*sinc(y1/D1).*exp(-(r1/(4*D1)).^2);
z = (1:n-1)*Dz/(n-1);
sg = exp(-(x1/(0.47*N*d1)).^16).*exp(-(y1/(0.47*N*d1)).^16);
t_vac = repmat(sg, [1 1 n]);
[xn yn Uvac] = ang_spec_multi_prop(pt, wvl, d1, dn, z, t_vac);
Uvac = Uvac.*exp(-1i*pi/(wvl*R)*(xn.^2+yn.^2));   % 准直
roi = hypot(xn/D2, yn/D2) < 2;
Ivac = abs(Uvac).^2;
fprintf('    ROI 辐照度 std/mean = %.2e；准直后相位 RMS = %.2e rad\n', ...
    std(Ivac(roi))/mean(Ivac(roi)), sqrt(mean(angle(Uvac(roi)).^2)));

% 自绘 fig-9-13（复现书图 9.7）
figure('Position', [50 50 950 400]);
subplot(1,2,1);
imagesc(2*xn(1,:)/D2, 2*yn(:,1)/D2, Ivac); axis image; colorbar;
caxis([0 max(Ivac(hypot(xn/D2,yn/D2)<2.5))]);
title('(a) 辐照度 [W/m^2]'); xlabel('2x_n/D_2'); ylabel('2y_n/D_2');
subplot(1,2,2);
imagesc(2*xn(1,:)/D2, 2*yn(:,1)/D2, angle(Uvac)); axis image; colorbar;
caxis([-pi pi]); colormap(gray);
title('(b) 相位 [rad]（准直后）'); xlabel('2x_n/D_2'); ylabel('2y_n/D_2');
sgt = suptitle('真空传播的模型点源（自绘，复现书图 9.7）');
saveas(gcf, fullfile(outdir, 'fig-9-13.png'));

%% ================= 验证 7：FT 屏 vs 次谐波屏（书 Listing 9.4 参数） =================
fprintf('===== 验证 7：FT 屏 vs 次谐波屏 =====\n');
Ns = 256; ds = 2/Ns; r0s = 0.1; L0s = 100; l0s = 0.01; nscr9 = 40;
D_ft = zeros(Ns); D_sh = zeros(Ns);
for ii = 1:nscr9
    ph_ft = ft_phase_screen(r0s, Ns, ds, L0s, l0s);
    [plo phi_] = ft_sh_phase_screen(r0s, Ns, ds, L0s, l0s);
    [d1m ~] = str_fcn2_ft(ph_ft, ds);   D_ft = D_ft + d1m;
    [d2m ~] = str_fcn2_ft(plo+phi_, ds); D_sh = D_sh + d2m;
end
D_ft = D_ft/nscr9;  D_sh = D_sh/nscr9;
% 径向平均（各向同性）
[nxg nyg] = meshgrid((0:Ns-1)*ds);
rg = hypot(nxg, nyg);
edges = linspace(0, 12*r0s, 121);
[~, ~, bin] = histcounts(rg(:), edges);
cprof = zeros(120,1); Dftp = cprof; Dshp = cprof;
for bb = 1:120
    m = (bin==bb);
    if any(m(:))
        cprof(bb) = mean(rg(m)); Dftp(bb) = mean(D_ft(m)); Dshp(bb) = mean(D_sh(m));
    end
end
keep = cprof > 0;
cprof = cprof(keep); Dftp = Dftp(keep); Dshp = Dshp(keep);
e_sh = max(abs(Dshp - 6.88*(cprof/r0s).^(5/3)))/max(6.88*(cprof/r0s).^(5/3));
e_ft = max(abs(Dftp - 6.88*(cprof/r0s).^(5/3)))/max(6.88*(cprof/r0s).^(5/3));
fprintf('    大间隔（8-12 r0）最大相对偏差：SH = %.1f%%，FT = %.1f%%\n', e_sh*100, e_ft*100);

figure('Position', [50 50 1050 420]);
subplot(1,2,1);
[plo phi_] = ft_sh_phase_screen(r0s, Ns, ds, L0s, l0s);
imagesc(-1, 1, -1, 1, plo+phi_); caxis([-20 20]); colormap(gray); colorbar;
title('(a) 典型次谐波相位屏 [rad]');
subplot(1,2,2);
plot(cprof/r0s, 6.88*(cprof/r0s).^(5/3), '-', 'Color', [.6 .6 .6], 'LineWidth', 1.5); hold on;
plot(cprof/r0s, Dshp, 'k-.', cprof/r0s, Dftp, 'k:');
xlabel('|\Delta r|/r_0'); ylabel('D_\phi(|\Delta r|) [rad^2]');
legend('理论（式 9.44）', 'Simulated SH', 'Simulated FT');
grid on; ylim([0 500]); xlim([0 12]);
saveas(gcf, fullfile(outdir, 'fig-9-11.png'));

%% ================= 验证 4-6：湍流仿真 + 验证（书 Listing 9.8/9.5.5） =================
fprintf('===== 验证 4：湍流仿真 40 次实现 =====\n');
l0 = 0; L0 = inf;
nreals = 40;
zt = [0 z];
alpha_t = zt/zt(end);
delta = (1-alpha_t)*d1 + alpha_t*dn;         % 各平面网格间距
MCF2 = zeros(N); Uout1 = []; phz_last = cell(1, nreals);
[xg yg] = meshgrid((-N/2:N/2-1)*dn);
mask = double(hypot(xg, yg) < D2/2);         % 望远镜口径掩模
for idxreal = 1:nreals
    t = zeros(N, N, n);
    for idxscr = 1:n
        [plo phi_] = ft_sh_phase_screen(r0scrn(idxscr), N, delta(idxscr), L0, l0);
        t(:,:,idxscr) = sg.*exp(1i*(plo+phi_));
        if idxscr == n; phz_last{idxreal} = plo+phi_; end
    end
    [xn yn Uout] = ang_spec_multi_prop(pt, wvl, d1, dn, z, t);
    Uout = Uout.*exp(-1i*pi/(wvl*R)*(xn.^2+yn.^2));
    MCF2 = MCF2 + corr2_ft(Uout, Uout, mask, dn);
    if idxreal == 1; Uout1 = Uout; end
end

% 自绘 fig-9-14（复现书图 9.8）
figure('Position', [50 50 950 400]);
subplot(1,2,1);
imagesc(2*xn(1,:)/D2, 2*yn(:,1)/D2, abs(Uout1).^2); axis image; colorbar;
caxis([0 prctile(abs(Uout1).^2(hypot(xn/D2,yn/D2)<1), 99)]); colormap(gray);
rectangle('Position', [-1 -1 2 2], 'Curvature', [1 1], 'EdgeColor', 'w');
title('(a) 辐照度（白圈 = 口径）'); xlabel('2x_n/D_2'); ylabel('2y_n/D_2');
subplot(1,2,2);
imagesc(2*xn(1,:)/D2, 2*yn(:,1)/D2, angle(Uout1)); axis image; colorbar;
caxis([-pi pi]); colormap(gray);
title('(b) 相位 [rad]（准直后）'); xlabel('2x_n/D_2'); ylabel('2y_n/D_2');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);
saveas(gcf, fullfile(outdir, 'fig-9-14.png'));

% 验证 5：末屏结构函数（书图 9.9）
fprintf('===== 验证 5：相位屏结构函数 =====\n');
r0_last = r0scrn(end);
Dacc = zeros(N);
for ii = 1:nreals
    [d2m mrN] = str_fcn2_ft(phz_last{ii}, delta(end));
    Dacc = Dacc + d2m;
end
Dacc = Dacc/nreals;
[nxg nyg] = meshgrid((0:N-1)*delta(end));
rg = hypot(nxg, nyg);
edges = linspace(0, 12*r0_last, 241);
[~, ~, bin] = histcounts(rg(:), edges);
cprof5 = zeros(240,1); Dprof = cprof5;
for bb = 1:240
    m = (bin==bb);
    if any(m(:)); cprof5(bb) = mean(rg(m)); Dprof(bb) = mean(Dacc(m)); end
end
keep = cprof5 > 0; cprof5 = cprof5(keep); Dprof = Dprof(keep);
err_scr = max(abs(Dprof - 6.88*(cprof5/r0_last).^(5/3)))/max(6.88*(cprof5/r0_last).^(5/3));
fprintf('    末屏（r0 = %.3f m）%d 块平均结构函数：最大相对偏差 %.1f%%\n', ...
    r0_last, nreals, err_scr*100);

% 验证 6：相干因子（书图 9.10）
fprintf('===== 验证 6：观测面相干因子 =====\n');
MCDOC2 = abs(MCF2)/abs(MCF2(N/2+1, N/2+1));
rg = hypot(xg, yg);
edges = linspace(0, 1.6*r0sw, 161);
[~, ~, bin] = histcounts(rg(:), edges);
cprof6 = zeros(160,1); mup = cprof6;
for bb = 1:160
    m = (bin==bb);
    if any(m(:)); cprof6(bb) = mean(rg(m)); mup(bb) = mean(MCDOC2(m)); end
end
keep = cprof6 > 0; cprof6 = cprof6(keep); mup = mup(keep);
mu_th = exp(-3.44*(cprof6/r0sw).^(5/3));
fprintf('    |mu_sim - mu_th| 平均 = %.3f（r0,sw = %.2f cm）\n', ...
    mean(abs(mup-mu_th)), r0sw*100);

figure('Position', [50 50 1050 420]);
subplot(1,2,1);
plot(cprof5/r0_last, 6.88*(cprof5/r0_last).^(5/3), '-', 'Color', [.6 .6 .6], 'LineWidth', 1.5); hold on;
plot(cprof5/r0_last, Dprof, 'k-.');
xlabel('|\Delta r|/r_0'); ylabel('D_\phi [rad^2]'); ylim([0 500]); xlim([0 12]);
legend('理论（式 9.44）', 'Simulated（40 块屏）'); grid on;
title('(a) 末平面屏幕平均结构函数（自绘，复现书图 9.9）');
subplot(1,2,2);
plot(cprof6/r0sw, mu_th, '-', 'Color', [.6 .6 .6], 'LineWidth', 1.5); hold on;
plot(cprof6/r0sw, mup, 'k.-');
xlabel('r_n/r_0'); ylabel('Coherence Factor'); ylim([0 1]); xlim([0 1.6]);
legend('theory', 'simulated'); grid on;
title('(b) 观测面相干因子（自绘，复现书图 9.10）');
saveas(gcf, fullfile(outdir, 'fig-9-15.png'));
fprintf('===== 工程 1 完成 =====\n');

%% ================= 局部函数 =================
function [xn yn Uout] = ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t)
% 书 Listing 9.1：分步束传播（多平面角谱传播 + 每平面相位屏 t）。
% t 为空（[]）时退化为真空传播。全局相位 e^{ikDz} 按书约定省略。
N = size(Uin, 1);
n1 = -N/2 : N/2 - 1;                       % MATLAB 1-based：零坐标在 N/2+1
[nx ny] = meshgrid(n1);
k = 2*pi/wvl;
nsq = nx.^2 + ny.^2;                       % 吸收边界（书 9-12 行，等价改写防溢出）
w = 0.47*N;
sg = exp(-(nsq/w^2).^8);
z = [0 z];
n = length(z);
Delta_z = z(2:n) - z(1:n-1);
alpha = z/z(n);
delta = (1-alpha)*delta1 + alpha*deltan;
m = delta(2:n)./delta(1:n-1);
x1 = nx*delta(1);  y1 = ny*delta(1);
r1sq = x1.^2 + y1.^2;
Q1 = exp(1i*k/2*(1-m(1))/Delta_z(1)*r1sq);
Uin = Uin.*Q1;
if ~isempty(t); Uin = Uin.*t(:,:,1); end
for idx = 1:n-1
    deltaf = 1/(N*delta(idx));
    fX = nx*deltaf;  fY = ny*deltaf;
    fsq = fX.^2 + fY.^2;
    Z = Delta_z(idx);
    Q2 = exp(-1i*pi^2*2*Z/m(idx)/k*fsq);
    Uin = ift2(Q2.*ft2(Uin/m(idx), delta(idx)), deltaf);
    Uin = Uin.*sg;
    if ~isempty(t); Uin = Uin.*t(:,:,idx+1); end
end
xn = nx*delta(n);  yn = ny*delta(n);
rnsq = xn.^2 + yn.^2;
Q3 = exp(1i*k/2*(m(n-1)-1)/(m(n-1)*Z)*rnsq);
Uout = Q3.*Uin;
end

function G = ft2(g, delta)
% 书 Listing 3.1：连续 FT 的 DFT 近似（居中输出）
N = size(g, 1);
G = fftshift(fft2(ifftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
% 书 Listing 3.2：连续 IFT 的 DFT 近似（居中输出）
N = size(G, 1);
g = fftshift(ifft2(ifftshift(G)))*(N*delta_f)^2;
end

function phz = ft_phase_screen(r0, N, delta, L0, l0)
% 书 Listing 9.2：FT 法相位屏（式 9.76-9.79）
del_f = 1/(N*delta);
fx = (-N/2:N/2-1)*del_f;
[fx fy] = meshgrid(fx);
f = hypot(fx, fy);
if l0 > 0; fm = 5.92/l0/(2*pi); else; fm = inf; end
if isfinite(L0); f0 = 1/L0; else; f0 = 0; end
PSD_phi = 0.023*r0^(-5/3)*exp(-(f/fm).^2)./(f.^2+f0^2).^(11/6);
PSD_phi(N/2+1, N/2+1) = 0;
cn = (randn(N)+1i*randn(N)).*sqrt(PSD_phi)*del_f;
phz = real(ift2(cn, 1));
end

function [phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0)
% 书 Listing 9.3：次谐波增广相位屏（式 9.81）
D = N*delta;
phz_hi = ft_phase_screen(r0, N, delta, L0, l0);
[x y] = meshgrid((-N/2:N/2-1)*delta);
phz_lo = zeros(size(phz_hi));
for p = 1:3
    del_f = 1/(3^p*D);
    fx = (-1:1)*del_f;
    [fx fy] = meshgrid(fx);
    f = hypot(fx, fy);
    if l0 > 0; fm = 5.92/l0/(2*pi); else; fm = inf; end
    if isfinite(L0); f0 = 1/L0; else; f0 = 0; end
    PSD_phi = 0.023*r0^(-5/3)*exp(-(f/fm).^2)./(f.^2+f0^2).^(11/6);
    PSD_phi(2,2) = 0;
    cn = (randn(3)+1i*randn(3)).*sqrt(PSD_phi)*del_f;
    SH = zeros(N);
    for ii = 1:9
        SH = SH + cn(ii)*exp(1i*2*pi*(fx(ii)*x+fy(ii)*y));
    end
    phz_lo = phz_lo + SH;
end
phz_lo = real(phz_lo) - mean(real(phz_lo(:)));
end

function [strfn mean_r] = str_fcn2_ft(phz, delta)
% 书 ch3 Listing 3.7 的等价重写（ch3 PDF 缺失）：无偏逐对差分结构函数。
% 注意不能用 2(var-C)：次谐波屏的大幅低频模使全局方差 != 逐对均方，
% 会在小间隔高估 D 数倍（Python 开发时实测）。
N = size(phz, 1);
q = zeros(2*N);  q(1:N, 1:N) = phz;
ones_w = zeros(2*N);  ones_w(1:N, 1:N) = 1;
Fp = fft2(q);  Fo = fft2(ones_w);  F2 = fft2(q.^2);
C = real(ifft2(Fp.*conj(Fp)));            % C[n]  = Σ_{m>=n} q[m]q[m-n]
S2b = real(ifft2(F2.*conj(Fo)));          % S2b[n] = Σ_{m>=n} q[m]^2
S2a = real(ifft2(Fo.*conj(F2)));          % S2a[n] = Σ_{m<N-n} q[m]^2
ov = (N-(0:N-1))' * (N-(0:N-1));  ov(ov==0) = 1;
strfn = (S2a(1:N,1:N) + S2b(1:N,1:N) - 2*C(1:N,1:N))./ov;
[X Y] = meshgrid((0:N-1)*delta);
mean_r = hypot(X, Y);
end

function [mcf mean_r] = corr2_ft(u1, u2, mask, delta)
% 书 ch3 Listing 3.10 的等价重写：FT 法互相关（相关定理）。
N = size(u1, 1);
U1 = ft2(u1.*mask, delta);
U2 = ft2(u2.*mask, delta);
mcf = ift2(U1.*conj(U2), 1/(N*delta));
[nxg nyg] = meshgrid((-N/2:N/2-1)*delta);
mean_r = hypot(nxg, nyg);
end
