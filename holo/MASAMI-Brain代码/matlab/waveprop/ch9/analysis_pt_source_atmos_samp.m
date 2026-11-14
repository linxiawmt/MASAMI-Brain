% ===== 书 Schmidt PM199 Listing 9.6 =====
% 例题（9.5.2 节）：带湍流的采样约束分析（式 9.86-9.90）。
% 前提：先运行 example_pt_source_atmos_setup.m（用到 D1、D2、wvl、Dz、
% r0sw、R 等变量）。c = 2 的 Coy 湍流展宽模型（式 9.84/9.85）：
%   D1' = D1 + c*wvl*Dz/r0,rev，D2' = D2 + c*wvl*Dz/r0；
% 本例 Cn2 沿路径常数，反向传播的 r0,rev = r0sw，故两处都用 r0sw。
% 约束 1/2 用 D1'、D2'；约束 3（数值判据）不受湍流影响、原样保留。
% 选点结果（原书图 9.6 白叉）：delta1 = deltan = 1 cm、N = 2^9 = 512、
% 最少 2 个平面（本例用 11 个平面/10 次部分传播来铺开湍流）。

c = 2;                               % 湍流敏感系数（2 -> ~97% 光，4 -> ~99%）
D1p = D1 + c*wvl*Dz/r0sw;            % 式 (9.84)：湍流模糊后的源尺寸
D2p = D2 + c*wvl*Dz/r0sw;            % 式 (9.85)：湍流模糊后的观测区尺寸

delta1 = linspace(0, 1.1*wvl*Dz/D2p, 100);
deltan = linspace(0, 1.1*wvl*Dz/D1p, 100);
% constraint 1（式 9.86：deltan 上界，直线）
deltan_max = -D2p/D1p*delta1 + wvl*Dz/D1p;
% constraint 3（式 9.88：与真空相同，斜率 1 + Dz/R）
d2min3 = (1+Dz/R)*delta1 - wvl*Dz/D1p;
d2max3 = (1+Dz/R)*delta1 + wvl*Dz/D1p;
[delta1 deltan] = meshgrid(delta1, deltan);
% constraint 2（式 9.87：N 下界，画 log2 N 等高线 -> 原书图 9.6）
N2 = (wvl * Dz + D1p*deltan + D2p*delta1) ...
    ./ (2 * delta1 .* deltan);

% constraint 4 / 部分传播检查：选中 delta1 = deltan = 1 cm、N = 512
d1 = 10e-3;
d2 = 10e-3;
N = 512;
d1*d2 * N / wvl      % 约束 4（式 7.59 的反面）：delta1*deltan*N >= wvl*Dz？
zmax = min([d1 d2])^2 * N / wvl   % 式 (9.89)：单次部分传播的最大距离
nmin = ceil(Dz / zmax) + 1        % 式 (9.90)：最少平面数（2；实际用 11）

% 绘图（对照原书图 9.6：log2 N 等高线 + 约束 1 虚线 + 约束 3 点划线 + 白叉）
figure;
[d1g d2g] = meshgrid(linspace(1e-4, 1.1*wvl*Dz/D2p, 200), ...
                     linspace(1e-4, 1.1*wvl*Dz/D1p, 200));
N2g = (wvl*Dz + D1p*d2g + D2p*d1g) ./ (2*d1g.*d2g);
contourf(d1g*1e3, d2g*1e3, log2(N2g), 4:13, 'LineStyle', 'none');
colorbar; hold on;
d1v = linspace(1e-4, 1.1*wvl*Dz/D2p, 200);
plot(d1v*1e3, (-D2p/D1p*d1v + wvl*Dz/D1p)*1e3, 'k--', 'LineWidth', 2);
plot(d1v*1e3, ((1+Dz/R)*d1v + wvl*Dz/D1p)*1e3, 'k-.', 'LineWidth', 2);
plot(10, 10, 'w*', 'MarkerSize', 12);
xlim([0 45]); ylim([0 65]); colorbar;
xlabel('\delta_1 [mm]'); ylabel('\delta_n [mm]');
title('约束 1, 2, 3（对照原书图 9.6，白叉 = 选中点 10 mm, 10 mm）');
