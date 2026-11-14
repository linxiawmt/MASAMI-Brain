% ===== 书 Schmidt PM199 Listing 9.8 =====
% 例题（9.5.4 节）：湍流仿真——每个实现重新抽 11 块相位屏（各屏 r0
% 来自 Listing 9.5 的约束优化解 r0scrn），传播、准直、用 corr2_ft 累积
% 观测面（口径内）的互相干函数；nreals 个实现平均后按中心值归一化得
% 相干因子 MCDOC2，与理论 exp[-3.44/2*(rn/r0,sw)^(5/3)]（式 9.32+9.44）
% 对比（原书图 9.10）。
% 注意：正文说"重复 40 次"，本清单印的是 nreals = 20（原书两处不一致，
% Python 版按 20 实现）；第 31 行续行符书印 4 个点"...."，应为"..."。
% 前提：Listings 9.5-9.7 已运行（r0scrn、pt、sg、z、delta1、deltan、n、
% N、wvl、R、D2、xn、yn 等变量在基底工作区）。
% 依赖：ft_sh_phase_screen.m、ang_spec_multi_prop.m、corr2_ft.m、circ.m。

l0 = 0;      % inner scale [m]
L0 = inf;    % outer scale [m]

zt = [0 z];                      % propagation plane locations
Delta_z = zt(2:n) - zt(1:n-1);   % propagation distances
% grid spacings（各平面网格间距，随 alpha 线性插值——生成相位屏要用
% 各屏所在平面的 delta，屏幕物理尺寸 N*delta(idx) 随之不同）
alpha = zt / zt(n);
delta = (1-alpha) * delta1 + alpha * deltan;

% initialize array for phase screens
phz = zeros(N, N, n);
nreals = 20;                     % number of random realizations
% initialize arrays for propagated fields,
% aperture mask, and MCF
Uout = zeros(N);
mask = circ(xn/D2, yn/D2, 1);    % 观测望远镜口径（直径 D2 的圆）
MCF2 = zeros(N);
sg = repmat(sg, [1 1 n]);
for idxreal = 1 : nreals         % loop over realizations
    idxreal
    % loop over screens
    for idxscr = 1 : 1 : n
        [phz_lo phz_hi] ...
            = ft_sh_phase_screen ...
            (r0scrn(idxscr), N, delta(idxscr), L0, l0);
        phz(:,:,idxscr) = phz_lo + phz_hi;
    end
    % simulate turbulent propagation（t = 吸收边界 .* exp(i*相位屏)）
    [xn yn Uout] = ang_spec_multi_prop(pt, wvl, ....
        delta1, deltan, z, sg.*exp(i*phz));      %（书印 "...."，应为 "..."）
    % collimate the beam
    Uout = Uout .* exp(-i*pi/(wvl*R)*(xn.^2+yn.^2));
    % accumulate realizations of the MCF
    MCF2 = MCF2 + corr2_ft(Uout, Uout, mask, deltan);
end

% modulus of the complex degree of coherence（按中心值归一化）
MCDOC2 = abs(MCF2) / (MCF2(N/2+1,N/2+1));

% 可视化（对照原书图 9.8/9.10）
figure;
subplot(1,2,1);
imagesc(xn(1,:)/D2, yn(:,1)/D2, abs(Uout).^2); axis image; colorbar;
xlabel('2x_n/D_2'); ylabel('2y_n/D_2'); title('湍流辐照度 [W/m^2]');
subplot(1,2,2);
imagesc(xn(1,:)/D2, yn(:,1)/D2, angle(Uout)); axis image; colorbar;
xlabel('2x_n/D_2'); ylabel('2y_n/D_2'); title('湍流相位（准直后）[rad]');

% 理论曲线（式 9.32 + 9.44：mu = exp[-D/2] = exp[-3.44/2*(r/r0,sw)^(5/3)]）
figure;
th = 0 : pi/180 : 2*pi;
plot(cos(th)*1.5, sin(th)*1.5, 'k:', 'LineWidth', 1); hold on;   % r/r0 = 1.5 参考圈
rad = sqrt(xn.^2 + yn.^2);
bins = linspace(0, 1.5*r0sw, 40);
[~, ~, ib] = histcounts(rad(:), bins);
mu_b = accumarray(ib, MCDOC2(:), [numel(bins)-1 1], @mean, NaN);
rc = (bins(1:end-1) + bins(2:end)) / 2;
mu_th = exp(-3.44/2*(rc/r0sw).^(5/3));
plot(rc/r0sw, mu_th, 'k-', 'LineWidth', 1.5);
plot(rc/r0sw, mu_b, 'k.');
xlabel('r_n / r_0'); ylabel('Coherence Factor');
legend('r/r_0 = 1.5', 'theory', 'simulated');
title('观测面相干因子（对照原书图 9.10）');
