%% ===== 书 Schmidt PM199 第 5 章 =====
%  成像系统与像差：Zernike 圆多项式与波前模式分解
%  （原书 5.1 节，Listings 5.1-5.3；与 python/waveprop/第5章_代码实战.ipynb
%    工程 1 数值行为一致，本机无 MATLAB 时以 Python 实跑结果为准）
%
%  原书要点：
%    广义光瞳函数 P = P_ap .* exp(i*2*pi*W)（式 5.1），波前差 W 以"波"为单位；
%    Noll 单指标 Zernike：Z_i = sqrt(2(n+1)) R_n^m(r) G_m(theta)（式 5.3/5.4），
%    径向 R_n^m 按 gamma 函数实现（式 5.5a，原书注释：gamma 比 factorial 快）；
%    单位圆上正交（式 5.7-5.9）-> 系数 = 正交投影/最小二乘（式 5.10-5.18），
%    W = Z*A（式 5.17），A = (Z'*Z)^{-1} Z'*W（式 5.18，MATLAB 里就是 Z\W）；
%    模式受限版 W' = sum_{i<=nZ} a_i Z_i（式 5.19-5.20）；
%    波前方差 sigma^2 = sum_{i>=2} a_i^2（式 5.21-5.26）。
%
%  原书勘误（本脚本按修正实现）：
%    (a) 式 (5.4)/(5.3) 的 m=0 情形印成 Z_i = R_n^0(r)——缺 sqrt(n+1)
%        归一化因子；Listing 5.1 第 9 行 sqrt(n+1)*zrf(n,0,r) 才是对的
%        （对照表 5.2 的 Z4 = sqrt(3)(2r^2-1) 可证）。
%    (b) 式 (5.7) 印成 1/(2n+1)——应为 1/(2(n+1))。n=0 即反例：
%        int_0^1 1*1*r dr = 1/2 而非 1；n=2 时数值积分 = 1/6 而非 1/5。
%    (c) 表 5.2 标题写"前 36 项"，实际列出 i = 1..37（含三级球差）。
%
%  验证（编号与 notebook 工程 1 一致）：
%   (0) 表 5.2 抽查：生成的 (i,n,m) 映射与书中多项式逐一对照；
%   (1) Listing 5.2 复现：W = 0.5*Z2 + 0.25*Z4 - 0.6*Z21，Z\W 精确恢复权重；
%   (2) Gram 矩阵正交性（孔径内求和）+ 式 (5.7) 勘误的数值积分复核 +
%       sigma^2 两种算法互证（式 5.21 离散 vs 式 5.26）；
%   (3) Listing 5.3 复现：随机相位屏 -> 200 阶投影 -> 模式受限 3/16/36/100 阶，
%       残差方差两种算法互证（差值 = 200 阶装不下的网格尺度内容）。
%
%  MATLAB/Python 差异（对应 notebook 内注释）：
%    - 线性索引：MATLAB 按列主序（原书 5.1.2 节点名 MATLAB 用 column-major），
%      logical 掩码索引天然给出列主序向量，与式 (5.14) 约定一致；numpy 布尔
%      掩码是行主序——顺序不同不影响最小二乘的解。
%    - Z\W 即式 (5.18) 的最小二乘解，等价 np.linalg.lstsq。
%    - 相位屏随机数：rng 与 numpy default_rng 实现不同，屏不逐点一致，
%      统计结论一致。
%    - 本地函数需 R2016b+；sum(...,'all') / max(...,[],'all') 需 R2018b+。
clear; clc; close all;

%% 验证 0：表 5.2 抽查
rr = linspace(0, 1, 41);
tt = linspace(0, 2*pi, 73);
[RR, TT] = meshgrid(rr, tt);
checks = {
    1,  @(r,t) ones(size(r)),                          'Z1  = 1 (piston)';
    2,  @(r,t) 2*r.*cos(t),                            'Z2  = 2r cos(theta)';
    4,  @(r,t) sqrt(3)*(2*r.^2-1),                     'Z4  = sqrt3(2r^2-1)';
    11, @(r,t) sqrt(5)*(6*r.^4-6*r.^2+1),              'Z11 = sqrt5(6r^4-6r^2+1)';
    22, @(r,t) sqrt(7)*(20*r.^6-30*r.^4+12*r.^2-1),    'Z22 = sqrt7(20r^6-...-1)';
    37, @(r,t) 3*(70*r.^8-140*r.^6+90*r.^4-20*r.^2+1), 'Z37 = 3(70r^8-...+1)'};
fprintf('===== 验证 0：表 5.2 抽查 =====\n');
for k = 1:size(checks, 1)
    err = max(abs(zernike(checks{k,1}, RR, TT) - checks{k,2}(RR, TT)), [], 'all');
    fprintf('    %-24s max|diff| = %.2e\n', checks{k,3}, err);
end

%% 验证 1：Listing 5.2 投影复现（模式 2/4/21）
N = 32;  L = 2;  delta = L / N;
[x, y] = meshgrid((-N/2 : N/2-1) * delta);
r = hypot(x, y);  theta = atan2(y, x);
ap = double(hypot(x, y) <= 1);            % circ(x,y,2)：直径 2 -> 单位圆孔径
idxAp = logical(ap);
z2  = zernike(2,  r, theta)  .* ap;
z4  = zernike(4,  r, theta)  .* ap;
z21 = zernike(21, r, theta)  .* ap;
W = 0.5*z2 + 0.25*z4 - 0.6*z21;           % 原书 Listing 5.2 第 16 行
Z = [z2(idxAp) z4(idxAp) z21(idxAp)];     % 式 (5.16)：列主序线性索引
A = Z \ W(idxAp);                         % 式 (5.18) 的 MATLAB 一行版
fprintf('===== 验证 1：Listing 5.2 投影复现 =====\n');
fprintf('    恢复系数 A = [%+.6f %+.6f %+.6f]（真值 [0.5 0.25 -0.6]）\n', A);
fprintf('    拟合残差 max = %.2e\n', max(abs(Z*A - W(idxAp)), [], 'all'));

%% 验证 2：正交性 + 式 (5.7) 勘误复核 + sigma^2 互证
nG = 12;
G = zeros(nG);
for ii = 1:nG
    Zi = zernike(ii, r, theta);
    for jj = ii:nG
        Zj = zernike(jj, r, theta);
        G(ii, jj) = sum(Zi(idxAp) .* Zj(idxAp)) * delta^2;
    end
end
G = G + G.' - diag(diag(G));
off = max(abs(G - diag(diag(G))), [], 'all') / pi;
diagErr = max(abs(diag(G)/pi - 1), [], 'all');
fprintf('===== 验证 2：正交性 + 式 (5.7) 勘误 =====\n');
fprintf('    N=32：Gram 对角/pi 最大偏差 %.2e，非对角最大 |G_ij|/pi %.2e\n', ...
    diagErr, off);
% 式 (5.7) 书印 1/(2n+1)：对 n=2 即 1/5；真实值 = 1/[2(n+1)] = 1/6
val = integral(@(q) (2*q.^2 - 1).^2 .* q, 0, 1);
fprintf('    int_0^1 (R_2^0)^2 r dr = %.6f = 1/6 = 1/[2(n+1)]\n', val);
fprintf('    书式 (5.7) 印的是 1/(2n+1) = %.6f —— 与表 5.2/Listing 5.1 矛盾，系印刷错误\n', 1/5);
% RMS 方差恒等式（式 5.21 vs 5.26）：W = 0.5 Z2 + 0.25 Z4 - 0.6 Z21
% 全部 i>=2，sigma^2 应 = 0.5^2 + 0.25^2 + 0.6^2 = 0.6725
Wsm = 0.5*z2 + 0.25*z4 - 0.6*z21;
sig2_direct = mean((Wsm(idxAp) - mean(Wsm(idxAp))).^2);   % 式 (5.21) 离散版
fprintf('    光滑波前 sigma^2：N=32 直接 = %.4f vs 系数平方和 = %.4f（相对差 %.1e）\n', ...
    sig2_direct, 0.6725, abs(sig2_direct - 0.6725)/0.6725);

%% 验证 3：Listing 5.3 模式分解/合成 + 残差方差恒等式
N3 = 128;  L3 = 2.0;  delta3 = L3/N3;  r0 = L3/20;
[x3, y3] = meshgrid((-N3/2 : N3/2-1) * delta3);
ap3 = hypot(x3, y3) <= 1;
idxA = logical(ap3);
phi = phase_screen_min(r0, N3, delta3, 0.6, 7);   % rad
W3 = phi / (2*pi) .* ap3;                  % rad -> 波（Listing 5.3 第 16 行）
r3 = hypot(x3, y3);  t3 = atan2(y3, x3);
Wvec = W3(idxA);
fprintf('===== 验证 3：Listing 5.3 模式分解/合成 =====\n');
fprintf('    相位屏：sigma_W = %.4f 波（孔径内 RMS；重定标目标 0.6 为全帧）\n', std(Wvec));
nModes = 200;
Zmat = zeros(numel(Wvec), nModes);
for q = 1:nModes
    tmp = zernike(q, r3, t3);
    Zmat(:, q) = tmp(idxA);
end
Amode = Zmat \ Wvec;                       % 式 (5.18)
fprintf('    200 阶投影：a1..a5 = [%+.4f, %+.4f, %+.4f, %+.4f, %+.4f]\n', Amode(1:5));
% QR 重正交版：Q 列正交，sum(Bq(nZ+1:end)^2) 精确等于残差均方
[Q, ~] = qr(Zmat, 0);
Bq = Q \ Wvec;
Nap3 = numel(Wvec);
nZ_list = [3 16 36 100];
resD = zeros(size(nZ_list));  resC = zeros(size(nZ_list));  resQ = zeros(size(nZ_list));
for q = 1:numel(nZ_list)
    nZ = nZ_list(q);
    resD(q) = mean((Wvec - Zmat(:,1:nZ)*Amode(1:nZ)).^2);
    resC(q) = sum(Amode(nZ+1:end).^2) / Nap3;
    resQ(q) = mean((Wvec - Q(:,1:nZ)*Bq(1:nZ)).^2);
    fprintf('    nZ=%3d：残差直接=%.6f vs 系数平方和=%.6f（相对差 %.1e）；QR 直接=%.6f\n', ...
        nZ, resD(q), resC(q), abs(resD(q)-resC(q))/resD(q), resQ(q));
end
fprintf('    （结论：sigma^2 = sum_{i>=2} a_i^2 在光滑波前上到离散正交精度；'...
    '相位屏上两曲线差恒定 = 200 阶装不下的网格尺度内容，对应式 (5.15) 的告诫）\n');

%% 绘图：模式受限屏（仿原书图 5.3/5.4）+ 系数谱 + 残差方差
scrData = {W3, '(a) 原始相位屏'};
lab = {'3', '16', '36', '100'};
tags = {'b', 'c', 'd', 'e'};
for q = 1:4
    scr = zeros(N3);
    scr(idxA) = Zmat(:, 1:nZ_list(q)) * Amode(1:nZ_list(q));
    scrData{end+1} = scr;                 %#ok<SAGROW>
    scrData{end+1} = sprintf('(%s) 限 %s 阶', tags{q}, lab{q});  %#ok<SAGROW>
end
figure('Position', [40 40 1400 640]);
for q = 1:5
    subplot(2, 5, q);
    imagesc(scrData{2*q-1} .* ap3);  colormap gray;  caxis([-1.2 1.2]);
    axis image off;  title(scrData{2*q});
end
subplot(2, 5, [6 7]);
semilogy(2:nModes, abs(Amode(2:end)), 'k.-', 'MarkerSize', 4);
set(gca, 'YTick', [1e0 1e-1 1e-2 1e-3 1e-4], 'YTickLabel', {'1','0.1','0.01','0.001','1e-4'});
xlabel('Zernike 阶号 i');  ylabel('|a_i|（波）');
title('(f) 相位屏的 Zernike 系数谱');  grid on;  ylim([1e-4 1]);
subplot(2, 5, [8 10]);
plot(nZ_list, resD, 'ko-', 'LineWidth', 1.2);  hold on;
plot(nZ_list, resC, 'r.--', 'LineWidth', 1.2, 'MarkerSize', 8);
plot(nZ_list, resQ, 'bx', 'MarkerSize', 11, 'LineStyle', 'none');
xlabel('保留模式数 nZ');  ylabel('残差方差 sigma^2（波^2）');
title('(g) 残差方差：两种算法互证');  grid on;
legend({'直接算 mean((W-W'')^2)', '系数平方和（式 5.26）', 'QR 重正交后直接算'}, ...
    'FontSize', 8, 'Location', 'northeast');

%% ===== 本地函数 =====
function [n, m] = zidx(i)
% 表 5.2 的 (i,n,m) 映射（Noll 单指标，i 从 1 开始；规则按表归纳）：
% 按 n 升序，第 n 行占 n+1 个模式；行内 |m| 与 n 同奇偶、升序；
% m=0 占 1 行，|m|>0 占 2 行（sin/cos 由 i 奇偶定：i 奇 sin，i 偶 cos，
% 即式 (5.5b)）。列主序/行主序无关紧要，映射本身逐项对上表 5.2。
persistent rows
if isempty(rows)
    rows = zeros(0, 2);
    for n0 = 0:40
        for m0 = mod(n0, 2):2:n0
            rows(end+1, :) = [n0, m0];    %#ok<AGROW>
            if m0 > 0
                rows(end+1, :) = [n0, m0];  %#ok<AGROW>
            end
        end
    end
end
n = rows(i, 1);
m = rows(i, 2);
end

function R = zrf(n, m, r)
% Zernike 径向多项式 R_n^m(r)（式 5.5a）。阶乘按原书注释用 gamma 实现
% （s! = Gamma(s+1)），比 factorial 快。
R = zeros(size(r));
for s = 0:(n-m)/2
    num = (-1)^s * gamma(n - s + 1);
    den = gamma(s+1) * gamma((n+m)/2 - s + 1) * gamma((n-m)/2 - s + 1);
    R = R + num / den * r.^(n - 2*s);
end
end

function Z = zernike(i, r, theta)
% 书式 (5.4)（m=0 分支已按勘误补 sqrt(n+1)）+ Listing 5.1。
% i 从 1 开始（i=1 是 piston）；r 是归一化半径（孔径边缘 r=1）。
[n, m] = zidx(i);
R = zrf(n, m, r);
if m == 0
    Z = sqrt(n+1) * R;                    % Listing 5.1 第 9 行
else
    if mod(i, 2) == 0                     % i 偶 -> cos（式 5.5b）
        Z = sqrt(2*(n+1)) * R .* cos(m*theta);
    else                                  % i 奇 -> sin
        Z = sqrt(2*(n+1)) * R .* sin(m*theta);
    end
end
end

function W = phase_screen_min(r0, N, delta, target_rms_waves, seed)
% 最小 Kolmogorov 谱法相位屏（单位 rad）。谱 Phi = 0.023 r0^{-5/3} f^{-11/3}
% （von Karman 去掉外尺度项），粗网格垫零 4 倍把最低频率降到 1/(4L) 弥补
% 无次谐波时的低频功率缺失，最后按孔径内 RMS 重定标到目标值（单位：波）。
% 正式的 ft_phase_screen（Kolmogorov + 次谐波）在第 9 章给出。
rng(seed, 'twister');
Npad = 4 * N;
fx = (0:Npad-1) / (Npad*delta);           % 与 ifft2 频率排列一致（= numpy fftfreq）
fx(fx >= 1/(2*delta)) = fx(fx >= 1/(2*delta)) - 1/delta;
[FX, FY] = meshgrid(fx, fx);
F2 = FX.^2 + FY.^2;
Phi = 0.023 * r0^(-5/3) .* F2.^(-11/6);
Phi(1, 1) = 0;                            % 活塞项（DC 处 0^{-11/6} 无定义）
W = real(ifft2(sqrt(Phi) .* fft2(randn(Npad))));
i0 = (Npad - N)/2 + 1;                    % 居中裁切（1-based）
W = W(i0:i0+N-1, i0:i0+N-1);
W = W / (std(W(:)) / (2*pi) / target_rms_waves);   % 重定标（单位：波）
end
