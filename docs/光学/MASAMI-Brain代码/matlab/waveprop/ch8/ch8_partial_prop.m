%% ===== 书 Schmidt PM199 第 8 章（工程 1）：用部分传播放宽采样约束 =====
%  原书 8.1-8.4 节（吸收边界、两次/任意次部分传播、采样约束），
%  书 Listing 8.1 / 8.2 原样整理 + 中文注释。
%  与 python/waveprop/ch8/ch8_partial_prop.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准。
%
%  原书要点：
%    角谱法只适合短距离：约束 4（N >= lambda*Dz/(delta1*delta2)，式 7.59）
%    随 Dz 线性变紧。解法：把全长传播拆成 n-1 次"部分传播"（n 个平面），
%    每段用角谱法，并在每个平面乘吸收边界压掉网格边缘能量、抑制 wrap-around：
%      - 吸收边界（8.1 节）：super-Gaussian（式 8.1）或 Tukey/余弦锥削窗
%        （式 8.2），中心 ~1、边缘 ~0；Hamming/Bartlett 中心也被压低，不适用。
%      - 网格间距沿 z 线性插值（式 8.8）：delta_i = (1-alpha_i)*delta1 + alpha_i*deltan，
%        alpha_i = z_i/Dz（相似三角形，图 8.3）。
%      - 中平面两个二次相位恰好互逆抵消（式 8.15），得干净算符链（式 8.16-8.18）。
%      - 采样约束（8.4 节）：约束 1/2 原样；约束 3 与分段数无关（式 8.22 ==
%        式 7.53，alpha 消去，习题 2）；约束 4 按段改写：
%        Delta_z_i <= min(delta1,deltan)^2*N/lambda（式 8.24），定出最少平面数。
%    书例（Listing 8.2，终稿参数）：方孔平面波 D1=2mm、lambda=1um、Dz=2m、
%    N=128、delta1=D1/30=66.7um、deltan=D2/30=200um、5 次部分传播（每段 0.4m）。
%
%  原书笔误/不一致（Python 侧数值 + 图面取证，详见笔记对应 warning）：
%    (a) 式 (8.2) Tukey 窗第一分支书印 "1, r >= alpha*L/2"：不等号印反
%        （平顶应在中心，与图 8.1 矛盾），且 r/L 与 alpha*N/2、N/2 单位混用；
%        按图 8.1 与标准 Tukey 定义实现（平顶 + 锥形）。
%    (b) 书例参数三处不一致：正文 delta_n<=133um（D2=4mm/30 来）vs
%        Listing 8.2 的 D2=6mm、deltan=200um；图 8.5 整张按旧稿 D2=4mm、
%        Dz=1m 绘制（约束 1 虚线截距 500um = lambda*1m/D1、斜率 -2 = -D2/D1，
%        白叉实测 (66.0,132.7)um 处旧稿口径 log2N=6.43 -> N>=87 -> 取 2 的幂
%        128，与正文"at least N=2^7=128"自洽）；图 8.6 标记间距 ~0.13mm 对应
%        133um。数值自洽的一套是 Listing 8.2（约束 2 需 N>=105；单次约束 4
%        需 N>=150，分段确实必要）；正文 133um 与 Dz=2m 组合违反约束 2
%        （N>=143）。以 Listing 8.2 为准复现。
%    (c) "至少 n = ceil(Dz/Dzmax)+1 次部分传播"的计数口径：n 按表 8.2 是
%        平面数（n-1 次传播），书例最少 4 次传播；Listing 8.2 取 5 次
%        （每段 0.4m），比下限多一段、更稳妥。
%    (d) 式 (8.16)（两次分段）只在中平面放 A[r2]，而式 (8.18)/Listing 8.1
%        在每个平面（含观察面）都乘吸收边界 —— 代码按式 (8.18) 实现。
%    (e) 表 8.2 把 alpha_i = z_i/Dz 描述成 "fractional distance from plane 1
%        to plane i+1"：按表 8.3 与 Listing 8.1，alpha_i 对应平面 i，措辞差一。
%
%  Python 实跑结果（本脚本按相同参数复现）：
%    验证 0：式 (8.15) 两中平面啁啾相对差 0.0；式 (8.20)/(8.22) 1000 组随机
%            参数判定 0 不一致（约束 3 与分段数无关）。
%    验证 1：n=1 单段（无吸收边界）与第 6 章 ang_spec_prop 最大复数差 0.0。
%    验证 2：Listing 8.2 复现，解析峰值 2.56（书图 8.6 峰值 2.57）、
%            第一暗环 1.00mm，ROI(+-3mm) 内 |dI|/pk = 3.9e-3、
%            相位残差 4.8e-3 rad；
%            正文 133um 口径 |dI|/pk = 4.3e-3（约束 2 突破仍大致可用）。
%    验证 4：Dz=8m 单次传播 ROI |dI|/pk = 1.12e-1，分段 20 次+吸收边界
%            8.2e-3 —— 差 14 倍（wrap-around 被抑制）。
%
%  MATLAB/Python 差异（对应 .py 内注释）：
%    - MATLAB 1-based（零坐标样本在 N/2+1 行）/ numpy 0-based（N//2 行）；
%    - MATLAB 全 double，无 int 溢出问题；numpy 侧 nsq.^8 必须转 float
%      （Windows int32 溢出坑，见 .py 的 super_gaussian_absorber 注释）；
%    - meshgrid/ft2/ift2 约定与书 Listing 2.5/2.6 一致。
clear; clc; close all;

%% ================= 验证 0：几何关系与 alpha 消去（式 8.8 / 8.15 / 习题 2） =====
fprintf('===== 验证 0：式 (8.8) 网格插值 / 式 (8.15) 相位抵消 / 式 (8.21)->(8.22) =====\n');
d1 = 66.7e-6;  d3 = 200e-6;  a1 = 0.3;
d2 = delta_middle(d1, d3, a1);                       % 式 (8.8)
fprintf('    式 (8.8)：delta1=66.7um, delta3=200um, alpha=0.3 -> delta2 = %.4f um\n', d2*1e6);
[lhs, rhs] = middle_chirp_cancels(d1, d3, 0.6, 1.4); % 式 (8.15)
fprintf('    式 (8.15)：(m2-1)/Dz2 = %.6e 与 (m1-1)/(m1*Dz1) = %.6e，相对差 %.1e\n', ...
    lhs, rhs, abs(lhs-rhs)/abs(rhs));
% 式 (8.21)->(8.22)：约束 3 的 alpha_2 消去（随机参数验证）
rng(8, 'twister');
mism = 0;
for it = 1:1000
    d1r = 10^(-5 + 2*rand);
    dnr = 10^(-5 + 2*rand);
    Dzr = 10^(-1 + 6*rand);
    Rr = 10^(-1 + 4*rand) * sign(rand - 0.5);
    wvlr = 10^(-7 + 2*rand);
    D1r = 10^(-3 + 2*rand);
    a2 = 0.05 + 0.9*rand;
    Dz1 = a2 * Dzr;
    d2r = (1 - a2)*d1r + a2*dnr;
    lo2 = (1 + Dz1/Rr)*d1r - wvlr*Dz1/D1r;           % 式 (8.20) 对 delta2
    hi2 = (1 + Dz1/Rr)*d1r + wvlr*Dz1/D1r;
    ok2 = (lo2 <= d2r) && (d2r <= hi2);
    lon = (1 + Dzr/Rr)*d1r - wvlr*Dzr/D1r;           % 式 (8.22) 对 deltan
    hin = (1 + Dzr/Rr)*d1r + wvlr*Dzr/D1r;
    okn = (lon <= dnr) && (dnr <= hin);
    if ok2 ~= okn
        mism = mism + 1;
    end
end
fprintf('    1000 组随机参数：式 (8.20)（对 delta2）与式 (8.22)（对 deltan）判定不一致 = %d 组\n', mism);

%% ================= 验证 2：Listing 8.2 书例复现（原书图 8.6） =================
fprintf('===== 验证 2：Listing 8.2 复现（方孔平面波，N=128，5 段传播） =====\n');
D1 = 2e-3;            % diameter of the source aperture [m]
D2 = 6e-3;            % diameter of the observation aperture [m]
wvl = 1e-6;           % optical wavelength [m]
z_tot = 2;            % propagation distance [m]
delta1 = D1/30;       % source-plane grid spacing [m]（30 点跨孔径）
deltan = D2/30;       % observation-plane grid spacing [m]（勘误 b：书正文 133um）
N = 128;              % number of grid points
nseg = 5;             % number of partial propagations
z = (1:nseg) * z_tot / nseg;      % switch from total distance to segments
[x1g y1g] = meshgrid((-N/2 : N/2-1) * delta1);
ap = rect(x1g/D1) .* rect(y1g/D1);                % source aperture
[xn yn Uout] = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z);
% analytic result for y=0 slice（MATLAB 1-based：零坐标在第 N/2+1 行）
Uan = fresnel_prop_square_ap(xn(N/2+1,:), 0, D1, wvl, z_tot);
win = abs(xn(N/2+1,:)) <= D2/2;
I_num = abs(Uout(N/2+1,:)).^2;  I_an = abs(Uan).^2;  pk = max(I_an);
eI = max(abs(I_num(win) - I_an(win)))/pk;
mph = win & (I_an > 0.05*pk);
eph = max(abs(angle(Uout(N/2+1,mph)./Uan(mph))));
fprintf('    delta1=%.2f um, deltan=%.2f um, N=%d, 段长 %.2f m x %d\n', ...
    delta1*1e6, deltan*1e6, N, z(1), nseg);
fprintf('    解析峰值辐照度 = %.2f（原书图 8.6 峰值 2.57）\n', pk);
fprintf('    ROI 内：|dI|/峰值 = %.2e；相位残差(|U|>5%%峰值) = %.2e rad\n', eI, eph);
% 正文口径（delta_n=133um）同样跑一遍：约束 2 被突破的后果
[xn1, ~, Uout133] = ang_spec_multi_prop_vac(ap, wvl, delta1, 133e-6, z);
Uan133 = fresnel_prop_square_ap(xn1(N/2+1,:), 0, D1, wvl, z_tot);
winb = abs(xn1(N/2+1,:)) <= 2e-3;
eI133 = max(abs(abs(Uout133(N/2+1,winb)).^2 - abs(Uan133(winb)).^2)) / max(abs(Uan133).^2);
mphb = winb & (abs(Uan133).^2 > 0.05*max(abs(Uan133).^2));
eph133 = max(abs(angle(Uout133(N/2+1,mphb)./Uan133(mphb))));
fprintf('    [正文口径 deltan=133um] +-2mm 内：|dI|/峰值 = %.2e；相位残差 = %.2e rad\n', eI133, eph133);

%% ================= 验证 3：采样约束总账 + 勘误 b 取证 =================
fprintf('===== 验证 3：书例采样约束 + 勘误 b（图 8.5 旧稿 / 正文 133um / Listing 200um）=====\n');
fprintf('    约束 1（D2=6mm, Dz=2m）：delta_n <= %.0f um（200um 满足）\n', ...
    constraint1_d2(delta1, D1, D2, wvl, z_tot)*1e6);
fprintf('    约束 2（deltan=200um, Dz=2m）：N >= %.1f = 2^%.2f -> N=128 成立\n', ...
    constraint2_N(delta1, 200e-6, D1, D2, wvl, z_tot), ...
    log2(constraint2_N(delta1, 200e-6, D1, D2, wvl, z_tot)));
fprintf('    约束 2（deltan=133um, D2=4mm, Dz=2m，正文口径）：N >= %.1f -> N=128 不成立\n', ...
    constraint2_N(delta1, 133e-6, D1, 4e-3, wvl, z_tot));
fprintf('    约束 2（deltan=133um, D2=4mm, Dz=1m，图 8.5 旧稿口径）：N >= %.1f = 2^%.2f -> N=128 成立\n', ...
    constraint2_N(delta1, 133e-6, D1, 4e-3, wvl, 1.0), ...
    log2(constraint2_N(delta1, 133e-6, D1, 4e-3, wvl, 1.0)));
[lo3, hi3] = constraint3_d2(delta1, D1, wvl, z_tot, inf);
fprintf('    约束 3（R=inf）：%.0f <= delta_n <= %.0f um（|1+Dz/R|=1 < D2/D1=3，免检）\n', ...
    lo3*1e6, hi3*1e6);
fprintf('    约束 4（单次，deltan=200um）：N >= %.1f > 128 -> 必须分段\n', ...
    wvl*z_tot/(delta1*deltan));
fprintf('    [图 8.5 取证] 白叉实测 (66.0, 132.7) um；虚线截距 500um=Dz=1m 旧稿、斜率 -2=D2=4mm 旧稿\n');
dzmax = dz_max_partial(delta1, deltan, N, wvl);   % 式 (8.24)/(8.25)
fprintf('    式 (8.25)：Delta_z_max = %.4f m（书印 0.567，尾数舍入）；最少平面数 n = ceil(2/%.4f)+1 = %d（%d 次传播）\n', ...
    dzmax, dzmax, ceil(z_tot/dzmax)+1, ceil(z_tot/dzmax));
fprintf('    Listing 8.2 取 5 次传播（每段 0.4 m）—— 比下限多一段，宁多勿少\n');

%% ================= 验证 4：Dz=8m wrap-around 对比（自绘 fig-8-7） =================
fprintf('===== 验证 4：长距传播的 wrap-around：单次 vs 分段 + 吸收边界 =====\n');
Dz4 = 8;
[~, ~, U_single] = ang_spec_prop(ap, wvl, delta1, deltan, Dz4);  % 单次、无吸收边界
Uan4 = fresnel_prop_square_ap(xn(N/2+1,:), 0, D1, wvl, Dz4);
I_an4 = abs(Uan4).^2;  pk4 = max(I_an4);
I_sg = abs(U_single(N/2+1,:)).^2;
nseg4 = 20;
z4 = (1:nseg4) * Dz4 / nseg4;
[~, ~, U_part4] = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z4);
I_pt4 = abs(U_part4(N/2+1,:)).^2;
eI_s4 = max(abs(I_sg(win) - I_an4(win)))/pk4;
eI_p4 = max(abs(I_pt4(win) - I_an4(win)))/pk4;
fprintf('    [Dz=8m] 约束 2 需 N >= %.0f、单次约束 4 需 N >= %.0f（均 > 128）\n', ...
    constraint2_N(delta1, deltan, D1, D2, wvl, Dz4), wvl*Dz4/(delta1*deltan));
fprintf('    单次：ROI |dI|/峰值 = %.2e；分段 20 次 + 吸收边界：%.2e —— 差 %.0f 倍\n', ...
    eI_s4, eI_p4, eI_s4/eI_p4);
ext_mm = 1e3 * deltan * N / 2;
figure('Position', [60 60 1350 440]);
subplot(1, 3, 1);
imagesc([-ext_mm ext_mm], [-ext_mm ext_mm], log10(abs(U_single) + 1e-12));
axis image; colorbar; colormap(inferno);
xlabel('x_n [mm]'); ylabel('y_n [mm]');
title('(a) 单次传播 Δz=8m（无吸收边界）');
subplot(1, 3, 2);
imagesc([-ext_mm ext_mm], [-ext_mm ext_mm], log10(abs(U_part4) + 1e-12));
axis image; colorbar; colormap(inferno);
xlabel('x_n [mm]'); ylabel('y_n [mm]');
title('(b) 分段 20 次 + 吸收边界');
subplot(1, 3, 3);
semilogy(xn(N/2+1,:), I_an4/pk4, 'k-', xn(N/2+1,:), I_pt4/pk4, 'r.', ...
    xn(N/2+1,:), I_sg/pk4, 'b.');
xlim([-3 3]); ylim([1e-5 2]); grid on;
xlabel('x_n [mm]'); ylabel('归一化辐照度');
legend('解析', '分段 20 次 + 吸收', '单次（无吸收）', 'Location', 'southwest');
title('(c) ROI 内辐照度切片（对数坐标）');

%% ================= 验证 5：各平面演化（自绘 fig-8-8） =================
fprintf('===== 验证 5：各平面场演化（吸收边界逐面吸收边缘能量） =====\n');
[~, ~, ~, planes] = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z, [], true);
nvec = -N/2 : N/2-1;
[nxs, nys] = meshgrid(nvec, nvec);
edge = sqrt(nxs.^2 + nys.^2) >= 0.45*N;
figure('Position', [60 60 1150 440]);
subplot(1, 2, 1);
hold on;
for i = 1:numel(planes)
    semilogy(nvec*delta1*1e3, abs(planes{i}(N/2+1,:)), 'LineWidth', 1.1, ...
        'DisplayName', sprintf('平面 %d（z=%.1f m）', i+1, z(i)));
end
plot(nvec*delta1*1e3, abs(ap(N/2+1,:)), 'k--', 'DisplayName', '源孔径');
set(gca, 'YScale', 'log'); ylim([1e-4 2]); grid on;
xlabel('索引 × δ_1 [mm]'); ylabel('|U|（y=0 切片）');
legend('Location', 'northeast');
title('(a) 各平面 |U| 截面：边缘被 super-Gaussian 逐面剪掉');
subplot(1, 2, 2);
ef = zeros(1, numel(planes));
for i = 1:numel(planes)
    ef(i) = sum(abs(planes{i}(edge)).^2)/sum(abs(planes{i}(~edge)).^2 ...
        + sum(abs(planes{i}(edge)).^2));
end
plot(2:numel(planes)+1, ef, 'ko-', 'MarkerFaceColor', 'k');
grid on;
xlabel('平面编号'); ylabel('边缘能量占比（r ≥ 0.45N/2）');
title('(b) 边缘能量占比被吸收边界逐面压制');

%% ================= 图 8-9：约束图复现（对照原书图 8.5） =================
d1g = linspace(1e-6, 260e-6, 300);
figure('Position', [80 80 820 560]);
[Nc2g, dng] = meshgrid(d1g, linspace(1e-6, 560e-6, 300));
C2 = constraint2_N(Nc2g, dng, D1, D2, wvl, z_tot);
contourf(Nc2g*1e6, dng*1e6, log2(max(C2, 1e-6)), 3:13, 'LineStyle', 'none');
colormap(gray); colorbar; hold on;
plot(d1g*1e6, constraint1_d2(d1g, D1, D2, wvl, z_tot)*1e6, 'k--', 'LineWidth', 2);
plot(d1g*1e6, constraint1_d2(d1g, D1, 4e-3, wvl, 1.0)*1e6, 'b--', 'LineWidth', 1.6);
plot(delta1*1e6, 200, 'r*', 'MarkerSize', 14);
plot(66.7, 133, 'w*', 'MarkerSize', 12);
xlim([0 260]); ylim([0 560]);
xlabel('\delta_1 [\mum]'); ylabel('\delta_n [\mum]');
title('图 8-9（自绘）约束图复现与勘误 b：图 8.5 按 D2=4mm、Δz=1m 旧稿绘制');
legend({'约束 2 等高线（log2 N，终稿）', '约束 1（终稿，截距 1000\mum)', ...
    '约束 1（图 8.5 实绘：截距 500\mum）', 'Listing 8.2：\delta_n=200\mum', ...
    '图 8.5 白叉实测 (66.0,132.7)\mum'});

%% ================= 本地函数 =================
function [xn, yn, Uout, planes] = ang_spec_multi_prop_vac(Uin, wvl, delta1, deltan, z, sg, record)
% 任意次部分传播（书 Listing 8.1 原样整理 + 可选输出）。
%   Uin    : U(r1)，源平面光场 [(W/m^2)^(1/2)]
%   wvl    : 光波长 [m]
%   delta1 : 源平面网格间距 [m]
%   deltan : 观察平面网格间距 [m]
%   z      : 平面位置数组 [z2, z3, ..., zn]（不含 z1=0）[m]
%   sg     : 吸收边界阵列（[] = 书默认 super-Gaussian w=0.47N、n=16）
%   record : true 时额外返回每个平面（乘吸收边界后）的场元胞数组
if nargin < 6 || isempty(sg)
    N0 = size(Uin, 1);
    sg = super_gaussian_absorber(N0);            % 书第 9-12 行
end
if nargin < 7
    record = false;
end
N = size(Uin, 1);                                % number of grid points
[nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1));       % 书：-N/2 : N/2-1（1-based）
k = 2*pi/wvl;                                    % optical wavevector
z = [0 z];                                       % propagation plane locations
n = length(z);
Delta_z = z(2:n) - z(1:n-1);                     % propagation distances
alpha = z / z(n);                                % fractional distances（alpha_i = z_i/Dz）
delta = (1-alpha) * delta1 + alpha * deltan;     % 式 (8.8)：各平面网格间距
m = delta(2:n) ./ delta(1:n-1);                  % scaling factors m_i
x1 = nx * delta(1);
y1 = ny * delta(1);
r1sq = x1.^2 + y1.^2;
Q1 = exp(i*k/2*(1-m(1))/Delta_z(1)*r1sq);        % 书第 26 行
Uin = Uin .* Q1;
planes = cell(1, n-1);
for idx = 1 : n-1                                % 书第 28-40 行
    % spatial frequencies (of i-th plane)
    deltaf = 1 / (N*delta(idx));
    fX = nx * deltaf;
    fY = ny * deltaf;
    fsq = fX.^2 + fY.^2;
    Z = Delta_z(idx);                            % propagation distance
    % quadratic phase factor（式 8.18 的 Q2[-Delta_z_i/m_i]，负号已对照原书）
    Q2 = exp(-i*pi^2*2*Z/m(idx)/k*fsq);
    % compute the propagated field（含 A[r_{i+1}] 吸收边界，按式 8.18）
    Uin = sg .* ift2(Q2 .* ft2(Uin / m(idx), delta(idx)), deltaf);
    if record
        planes{idx} = Uin;
    end
end
% observation-plane coordinates
xn = nx * delta(n);
yn = ny * delta(n);
rnsq = xn.^2 + yn.^2;
% 循环结束后 Z 恰为最后一段 Delta_z(n-1)（书代码同款写法）
Q3 = exp(i*k/2*(m(n-1)-1)/(m(n-1)*Z)*rnsq);
Uout = Q3 .* Uin;
end

function g = super_gaussian_absorber(N, w_frac, order)
% super-Gaussian 吸收边界（式 8.1；书 Listing 8.1 第 10-12 行的函数化）。
% sg = exp(-nsq.^8 / w^16)，nsq = nx^2+ny^2（索引单位）、w = 0.47*N，
% 即 n=16、sigma = 0.47*N*delta = 0.47*L（L 为网格全宽）。
if nargin < 2
    w_frac = 0.47;
end
if nargin < 3
    order = 16;
end
[nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1));
nsq = nx.^2 + ny.^2;
w = w_frac * N;
g = exp(-nsq.^(order/2) / w^order);
end

function g = tukey_absorber(N, alpha)
% Tukey（余弦锥削）吸收边界，式 (8.2) 的修正实现（勘误 a）。
% 书印第一分支 "r >= alpha*L/2" 不等号反了（平顶应在中心），且 r 与 r/L、
% L 与 N 单位混用；按图 8.1 与标准 Tukey 定义：平顶 r <= alpha*L/2，
% 锥形 alpha*L/2 <= r <= L/2（L = N*delta 全宽，r 取索引半径 * delta）。
% 方格四角（索引半径 > N/2）超出书公式定义域，钳为 0。
[nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1));
R = sqrt(nx.^2 + ny.^2);                         % 索引半径（单位：格点）
half = N/2;
a_half = alpha * half;
g = ones(N);
band = (R > a_half) & (R <= half);
g(band) = 0.5 * (1 + cos(pi*(R(band) - a_half)/((1-alpha)*half)));
g(R > half) = 0;
end

function d2 = delta_middle(d1, d3, alpha)
% 中间平面网格间距（式 8.8）：delta2 = (1-alpha)*delta1 + alpha*delta3。
d2 = (1 - alpha)*d1 + alpha*d3;
end

function [lhs, rhs] = middle_chirp_cancels(d1, d3, Dz1, Dz2)
% 式 (8.15)：(m2-1)/Dz2 == (m1-1)/(m1*Dz1)（中平面两个 Q 因子互逆）。
a = Dz1/(Dz1 + Dz2);
d2 = delta_middle(d1, d3, a);
m1 = d2/d1;
m2 = d3/d2;
lhs = (m2 - 1)/Dz2;
rhs = (m1 - 1)/(m1*Dz1);
end

function out = constraint1_d2(d1, D1, D2, wvl, Dz)
% 约束 1（式 7.14）：delta_n 上界（几何，不随分段改变）。
out = -D2/D1.*d1 + wvl*Dz/D1;
end

function Nlo = constraint2_N(d1, dn, D1, D2, wvl, Dz)
% 约束 2（式 7.20）：网格点数下界（几何，不随分段改变）。
Nlo = D1./(2*d1) + D2./(2*dn) + wvl*Dz./(2*d1.*dn);
end

function [lo, hi] = constraint3_d2(d1, D1, wvl, Dz, R)
% 约束 3（式 8.22 == 式 7.53）：与分段数无关（alpha 消去，习题 2）。
slope = 1 + Dz/R;
off = wvl*Dz/D1;
lo = slope*d1 - off;
hi = slope*d1 + off;
end

function dzmax = dz_max_partial(d1, dn, N, wvl)
% 式 (8.24)：单段最大传播距离 Delta_z_max = min(d1,dn)^2 * N / lambda。
dzmax = min(d1, dn)^2 * N / wvl;
end

function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
% 角谱传播（书 Listing 6.5；第 6 章），此处用于"单次传播"对照。
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
r1sq = x1.^2 + y1.^2;
df1 = 1/(N*d1);
[fX, fY] = meshgrid(n*df1, n*df1);
fsq = fX.^2 + fY.^2;
m = d2/d1;
[x2, y2] = meshgrid(n*d2, n*d2);
r2sq = x2.^2 + y2.^2;
Q1 = exp(1i*kk/2*(1 - m)/Dz*r1sq);
Q2 = exp(-1i*pi^2*2*Dz/m/kk*fsq);
Q3 = exp(1i*kk/2*(m - 1)/(m*Dz)*r2sq);
Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin/m, d1), df1);
end

function U = fresnel_prop_square_ap(x2, y2, D, wvl, Dz)
% 边长 D 方孔的 Fresnel 衍射解析解（书式 1.60/附录 B.5 思路；省略 e^{ikDz}）。
s = sqrt(2/(wvl*Dz));
[tX1, tX2] = deal(s*(x2 - D/2), s*(x2 + D/2));
[tY1, tY2] = deal(s*(y2 - D/2), s*(y2 + D/2));
[Sx1, Cx1] = fresnel_int(tX1);
[Sx2, Cx2] = fresnel_int(tX2);
[Sy1, Cy1] = fresnel_int(tY1);
[Sy2, Cy2] = fresnel_int(tY2);
Ix = sqrt(wvl*Dz/2)*((Cx2 - Cx1) + 1i*(Sx2 - Sx1));
Iy = sqrt(wvl*Dz/2)*((Cy2 - Cy1) + 1i*(Sy2 - Sy1));
U = Ix.*Iy/(1i*wvl*Dz);
end

function [S, C] = fresnel_int(t)
% Fresnel 积分 S(t)/C(t)（MATLAB 基础版无内置，用 integral() 逐点求值；
% Python 版用 scipy.special.fresnel）。
S = arrayfun(@(tt) integral(@(u) sin(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
C = arrayfun(@(tt) integral(@(u) cos(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）。
G = fftshift(fft2(fftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
% 原书 Listing 2.6：定标正确的二维逆 DFT（乘 (N*delta_f)^2）。
N = size(G, 1);
g = fftshift(ifft2(fftshift(G)))*(N*delta_f)^2;
end

function out = rect(x)
% rect(x)：|x|<1/2 取 1，|x|=1/2 取 1/2，其余 0。
out = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end
