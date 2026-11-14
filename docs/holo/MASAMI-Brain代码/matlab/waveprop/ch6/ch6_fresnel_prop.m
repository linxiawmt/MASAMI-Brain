%% ===== 书 Schmidt PM199 第 6 章（工程 1）=====
%  真空中的菲涅耳衍射：单步 / 两步 / 角谱 / ABCD 四种传播器
%  （原书 6.1-6.5 节，Listings 6.1-6.7）
%  与 python/waveprop/ch6/ch6_fresnel_prop.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准。
%
%  原书要点：
%    Fresnel 积分（式 6.1）的两种数值形态：
%      单步 FT 形态（式 6.5，Listing 6.1）：U = Q V F Q —— 一次 FT，
%        观测面采样 delta2 = lambda*Dz/(N*delta1) 固定（式 6.16）；
%      卷积形态（式 6.6）：两步传播（式 6.17-6.30，Listing 6.3）与
%        角谱传播（式 6.31/6.50，Listing 6.5）都引入缩放参数 m = delta2/delta1。
%    算子记号（式 6.7-6.13）：Q/V/F/F^-1/R/Q2，贯穿第 6-8 章。
%    简单光学系统（式 6.68-6.81）：ABCD 阵 + 广义 Huygens-Fresnel 积分
%    （Collins，式 6.77）的 Lambert-Fraser 卷积化（式 6.79-6.81）；
%    Listing 6.6 实现其自由 m 版本（笔记补出算子链，见验证 5）。
%
%  原书勘误（本脚本按修正实现并数值验证，详见笔记）：
%    (a) 式 (6.66) 的算子链不是 (6.50) 的正确 m->-m 代换：Q 系数应为
%        (m+1)/(m*Dz) 与 (1+m)/Dz，F^-1 求值点应为 -r2/m（(6.67) 顶符号才对）；
%        实输入场下书印链输出 = -conj(真值)（Python 实测 4e-16）。
%    (b) Listing 6.3 对负距离步用 |Dz1| 建中间面坐标网格，坐标对应翻转：
%        m != 1 时 "-" 支输出是真值的镜像（对称测试场不可见）。
%
%  验证（编号与 Python 版一致）：
%   (0) ft2/ift2 往返 + Parseval；
%   (1) Listing 6.2（N=1024 方孔 D=2mm, Dz=1m）：one_step vs 解析闭式；
%   (2) Listing 6.4（m=10.24）：two_step vs one_step vs 解析；表 6.2、
%       式 (6.19)-(6.23)、(6.55)/(6.56) 核对；
%   (3) 非对称场（m=2.5）：镜像效应实锤 + 修正版 == ang_spec（1e-14）；
%   (4) 式 (6.66) 勘误：书印链 = -conj(真值)；
%   (5) Listing 6.7：ABCD（f=inf）== ang_spec；ABCD（薄透镜 f=2m）==
%       ang_spec(孔径 x 发散相位)（同一离散算子）；
%   (6) 二次相位采样判据。
%
%  MATLAB/Python 差异（对应 .py 内注释）：
%    - MATLAB 1-based / numpy 0-based：居中排列的零坐标样本位置不同；
%    - meshgrid 行为一致；cart2pol -> hypot/atan2；
%    - 偶数网格的居中坐标反转 P(-r)：circshift(rot90(P,2),[1 1])
%      （numpy: roll(flip(P),(1,1))）；
%    - MATLAB 无内置 Fresnel 积分（符号工具箱除外），本脚本用 integral()
%      数值积分实现（Python 版用 scipy.special.fresnel）；
%    - 本地函数需 R2016b+；sum(...,'all') / max(...,[],'all') 需 R2018b+。
clear; clc; close all;

%% 验证 0：ft2/ift2 往返 + Parseval
N0 = 64;  L0 = 1.0;  d0 = L0/N0;  df0 = 1/(N0*d0);
rng(0, 'twister');
g0 = randn(N0) + 1i*randn(N0);
rt = ift2(ft2(g0, d0), df0);
G0 = ft2(g0, d0);
Eg = sum(abs(g0(:)).^2)*d0^2;
EG = sum(abs(G0(:)).^2)*df0^2;
fprintf('===== 验证 0：ft2/ift2 往返 + Parseval =====\n');
fprintf('    max|ift2(ft2(g)) - g| = %.2e\n', max(abs(rt(:) - g0(:))));
fprintf('    Parseval 相对差 = %.1e\n', abs(Eg - EG)/Eg);

%% 原书 Listing 6.2/6.4 的参数（方形孔径，Dz = 1 m，N = 1024）
N  = 1024;      % number of grid points per side
L  = 1e-2;      % total size of the grid [m]
d1 = L/N;       % grid spacing [m]
D  = 2e-3;      % diameter of the aperture [m]
wvl = 1e-6;     % optical wavelength [m]
k  = 2*pi/wvl;
Dz = 1;         % propagation distance [m]
n  = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
ap = rect(x1/D) .* rect(y1/D);

%% 验证 1：Listing 6.2 单步传播 vs 解析（原书图 6.1）
[x2, y2, U1] = one_step_prop(ap, wvl, d1, Dz);
Uan = fresnel_prop_square_ap(x2, y2, D, wvl, Dz);
win = (abs(x2) <= 3e-3) & (abs(y2) <= 3e-3);
I1  = abs(U1).^2;  Ian = abs(Uan).^2;
pkA = max(Ian(win));
ratio = U1(win)./Uan(win);
strong = abs(Uan(win)) > 0.05*pkA;
fprintf('===== 验证 1：Listing 6.2 单步传播 vs 解析 =====\n');
fprintf('    delta2 = %.4f mm（式 6.16），观测网格半宽 %.1f mm\n', ...
    wvl*Dz/(N*d1)*1e3, max(x2(:))*1e3);
fprintf('    轴上峰值辐照度：数值 %.4f vs 解析 %.4f（相对差 %.1e）\n', ...
    I1(N/2+1, N/2+1), Ian(N/2+1, N/2+1), ...
    abs(I1(N/2+1, N/2+1) - Ian(N/2+1, N/2+1))/Ian(N/2+1, N/2+1));
fprintf('    ±3 mm 窗内：max|U1-Uan|/max|Uan| = %.2e\n', ...
    max(abs(U1(win) - Uan(win)))/max(abs(Uan(win))));
fprintf('    |U|>5%%峰值处：|比值-1| max = %.2e；相位残差 max = %.2e rad\n', ...
    max(abs(abs(ratio(strong)) - 1)), ...
    max(abs(angle(ratio(strong)) - mean(angle(ratio(strong))))));

%% 验证 2：Listing 6.4 两步传播（m=10.24）+ 表 6.2 + 网格关系
d2_book = wvl*Dz/(N*d1);          % 原书取 one-step 的 delta2 -> m = 10.24
m_book = d2_book/d1;
[~, ~, U2] = two_step_prop(ap, wvl, d1, d2_book, Dz);
[~, ~, U2p] = two_step_prop_plus(ap, wvl, d1, d2_book, Dz);
I2 = abs(U2).^2;
fprintf('===== 验证 2：Listing 6.4 两步传播（m=10.24） =====\n');
fprintf('    m = %.2f：Dz1 = Dz/(1-m) = %.4f m（<0，图 6.2 几何）\n', ...
    m_book, Dz/(1 - m_book));
fprintf('    轴上峰值辐照度：two_step %.4f vs 解析 %.4f（相对差 %.1e）\n', ...
    I2(N/2+1, N/2+1), Ian(N/2+1, N/2+1), ...
    abs(I2(N/2+1, N/2+1) - Ian(N/2+1, N/2+1))/Ian(N/2+1, N/2+1));
fprintf('    书印"-"支 vs one_step：%.2e（方孔对称，镜像不可见）；"+"支 vs one_step：%.2e\n', ...
    max(abs(U2(:) - U1(:)))/max(abs(U1(:))), ...
    max(abs(U2p(:) - U1(:)))/max(abs(U1(:))));
m_t = 2.0;  d2_t = m_t*d1;
Dz1m = Dz/(1 - m_t);  Dz1p = Dz/(1 + m_t);
d1a_m = wvl*abs(Dz1m)/(N*d1);
fprintf('    m=2：表 6.2：Dz1-/Dz = %.1f（书 -1），Dz2-/Dz = %.1f（书 2）；', ...
    Dz1m/Dz, (Dz - Dz1m)/Dz);
fprintf('Dz1+/Dz = %.4f（书 1/3），Dz2+/Dz = %.4f（书 2/3）\n', ...
    Dz1p/Dz, (Dz - Dz1p)/Dz);
fprintf('    式 (6.19)-(6.23)：delta1a(-) = %.2f um；delta2 复得 = %.2f um = m*delta1\n', ...
    d1a_m*1e6, wvl*abs(Dz - Dz1m)/(N*d1a_m)*1e6);
fprintf('    式 (6.55)/(6.56)：1/(1-m)=%.6f vs d1/(d1-d2)=%.6f；m/(1-m)=%.6f vs d2/(d1-d2)=%.6f\n', ...
    1/(1 - m_t), d1/(d1 - d2_t), m_t/(1 - m_t), d2_t/(d1 - d2_t));

%% 验证 3：非对称场（m=2.5）：镜像效应与离散路径
Na = 512;  L1a = 8e-3;  d1a_ = L1a/Na;
m_a = 2.5;  d2a_ = m_a*d1a_;  Dza = 0.7;
na = -Na/2 : Na/2-1;
[xa, ya] = meshgrid(na*d1a_, na*d1a_);
ap_asym = rect((xa - 0.8e-3)/1.2e-3).*rect(ya/2e-3) ...
    + 0.6*circ_(xa - 1.8e-3, ya + 1.2e-3, 0.5e-3);
[~, ~, U_asym] = ang_spec_prop(ap_asym, wvl, d1a_, d2a_, Dza);
[~, ~, U_mc] = ang_spec_mbranch(ap_asym, wvl, d1a_, d2a_, Dza, true);
[~, ~, U_mb] = ang_spec_mbranch(ap_asym, wvl, d1a_, d2a_, Dza, false);
[~, ~, U_tm] = two_step_prop(ap_asym, wvl, d1a_, d2a_, Dza);
[~, ~, U_tp] = two_step_prop_plus(ap_asym, wvl, d1a_, d2a_, Dza);
[~, ~, U_tf] = two_step_prop_fixed(ap_asym, wvl, d1a_, d2a_, Dza);
pkA2 = max(abs(U_asym(:)));
U_asym_flip = rollflip(U_asym);
U_tm_flip = rollflip(U_tm);
fprintf('===== 验证 3：非对称场（m=2.5） =====\n');
fprintf('    [镜像勘误 b] 书印"-"支 vs ang_spec：%.2e；vs 镜像(ang_spec) |U| 差：%.2e\n', ...
    max(abs(U_tm(:) - U_asym(:)))/pkA2, ...
    max(abs(abs(U_tm(:)) - abs(U_asym_flip(:))))/pkA2);
fprintf('    修正"-"支（翻转中间场）vs ang_spec：%.2e —— 逐位一致\n', ...
    max(abs(U_tf(:) - U_asym(:)))/pkA2);
fprintf('    [离散路径] ang_spec vs (6.66)修正链：%.2e；"+"支 vs ang_spec：%.2e\n', ...
    max(abs(U_asym(:) - U_mc(:)))/pkA2, max(abs(U_tp(:) - U_asym(:)))/pkA2);
for m_scan = [1.5 2.0 2.5 3.0 4.0]
    [~, ~, A_] = ang_spec_prop(ap_asym, wvl, d1a_, m_scan*d1a_, Dza);
    [~, ~, C_] = ang_spec_mbranch(ap_asym, wvl, d1a_, m_scan*d1a_, Dza, true);
    fprintf('        m=%.1f：|mcorr-ang|/pk = %.3f\n', m_scan, ...
        max(abs(C_(:) - A_(:)))/max(abs(A_(:))));
end
m0 = (wvl*Dza/(Na*d1a_))/d1a_;
[~, ~, U_one_m0] = one_step_prop(ap_asym, wvl, d1a_, Dza);
[~, ~, U_tp_m0] = two_step_prop_plus(ap_asym, wvl, d1a_, m0*d1a_, Dza);
[~, ~, U_tf_m0] = two_step_prop_fixed(ap_asym, wvl, d1a_, m0*d1a_, Dza);
fprintf('    同网格互证（m0=%.2f）：one_step vs "+"支 = %.2e；vs 修正"-"支 = %.2e\n', ...
    m0, max(abs(U_one_m0(:) - U_tp_m0(:)))/max(abs(U_tp_m0(:))), ...
    max(abs(U_one_m0(:) - U_tf_m0(:)))/max(abs(U_tf_m0(:))));

%% 验证 4：式 (6.66) 勘误（书印链 = -conj(真值)）
fprintf('===== 验证 4：式 (6.66) 勘误 =====\n');
fprintf('    修正链 vs ang_spec：%.2e；书印链 vs ang_spec：%.2e\n', ...
    max(abs(U_mc(:) - U_asym(:)))/pkA2, max(abs(U_mb(:) - U_asym(:)))/pkA2);
fprintf('    书印链 vs -conj(ang_spec)：max|diff|/max|U| = %.2e（实输入场逐点恒等）\n', ...
    max(abs(U_mb(:) + conj(U_asym(:))))/pkA2);

%% 验证 5：Listing 6.7 的 ABCD 传播
ABCD_free = [1 Dz; 0 1];
[~, ~, U_af] = ang_spec_propABCD(ap, wvl, d1, d2_book, ABCD_free);
[~, ~, U_as] = ang_spec_prop(ap, wvl, d1, d2_book, Dz);
fprintf('===== 验证 5：Listing 6.7 的 ABCD 传播 =====\n');
fprintf('    f=inf（自由空间 ABCD）vs ang_spec：max|diff|/max|U| = %.2e\n', ...
    max(abs(U_af(:) - U_as(:)))/max(abs(U_as(:))));
f_lens = 2.0;
ABCD_lens = [1 Dz; 0 1]*[1 0; -1/f_lens 1];
d2_lens = 1.5*d1;
[~, ~, U_al] = ang_spec_propABCD(ap, wvl, d1, d2_lens, ABCD_lens);
[~, ~, U_lr2] = ang_spec_prop(ap.*exp(-1i*k*(x1.^2 + y1.^2)/(2*f_lens)), ...
    wvl, d1, d2_lens, Dz);
fprintf('    f=2 m 薄透镜 + Dz=1 m：ABCD 链 vs ang_spec(孔径 x e^{-ik r1^2/2f})：%.2e（同一离散算子）\n', ...
    max(abs(U_al(:) - U_lr2(:)))/max(abs(U_lr2(:))));

%% 验证 6：二次相位采样判据
fprintf('===== 验证 6：二次相位采样判据 =====\n');
fprintf('    源面啁啾（边缘 r=%.1f mm）：局部频率 %.0f /m vs Nyquist %.0f /m（比值 %.3f）\n', ...
    L/2*1e3, (L/2)/(wvl*Dz), 1/(2*d1), (L/2)/(wvl*Dz)*2*d1);
fprintf('    缩放链输入啁啾判据 |1-m|*r_edge/(lambda Dz)/(1/(2 d1))：m=2.5 时 %.3f\n', ...
    abs(1 - m_a)*(L1a/2)/(wvl*Dza)*2*d1a_);
fprintf('    (6.66)修正链曲率 (1+m)：比值 %.3f —— 更陡，是链间差的主因\n', ...
    (1 + m_a)*(L1a/2)/(wvl*Dza)*2*d1a_);

%% 绘图（对应 Python 版图 6-13）
sl = N/2 + 1;
m_ = abs(x2(1, :)) <= 3e-3;
figure('Position', [40 40 1500 820]);
subplot(2, 3, 1);
imagesc(x1(1, :)*1e3, y1(:, 1)*1e3, abs(ap)); axis image off;
title('(a) 源面：2 mm 方孔（Listing 6.2，N=1024）');
subplot(2, 3, 2);
plot(x2(1, m_)*1e3, Ian(sl, m_), 'k-', x2(1, m_)*1e3, I1(sl, m_), 'r.', ...
    x2(1, m_)*1e3, I2(sl, m_), 'b.');
title('(b) y2=0 辐照度切片：数值 vs 解析'); grid on;
xlabel('x2 [mm]'); ylabel('辐照度 [a.u.]');
subplot(2, 3, 3);
slw = m_ & (Ian(sl, :) > 0.05*max(Ian(sl, m_)));
ph1_ = angle(U1(sl, slw)./Uan(sl, slw));
ph1_ = ph1_ - mean(ph1_);
plot(x2(1, slw)*1e3, ph1_, 'b.');
title(sprintf('(c) one_step/解析 相位残差 max %.1e rad', max(abs(ph1_))));
xlabel('x2 [mm]'); ylabel('相位差 [rad]'); grid on;
subplot(2, 3, 4);
imagesc(xa(1, :)*1e3, ya(:, 1)*1e3, abs(U_asym)); axis image off;
title('(d) 非对称场 ang_spec（m=2.5, Dz=0.7 m）');
subplot(2, 3, 5);
sa = Na/2 + 1;
plot(xa(1, :)*1e3, abs(U_asym(sa, :))/pkA2, 'k-', ...
    xa(1, :)*1e3, abs(U_tm(sa, :))/pkA2, 'b--', ...
    xa(1, :)*1e3, abs(U_tm_flip(sa, :))/pkA2, 'r.');
title('(e) 勘误 b：书印两步"-"支输出镜像'); grid on;
xlabel('x2 [mm]'); ylabel('归一化 |U|');
subplot(2, 3, 6);
plot(xa(1, :)*1e3, real(U_asym(sa, :))/pkA2, 'k-', ...
    xa(1, :)*1e3, -real(conj(U_asym(sa, :)))/pkA2, 'b--', ...
    xa(1, :)*1e3, real(U_mc(sa, :))/pkA2, 'r.');
title('(f) 勘误 a：书印 (6.66) 链 = -conj(真值)'); grid on;
xlabel('x2 [mm]'); ylabel('归一化 Re U');

%% ===== 本地函数 =====
function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）
G = fftshift(fft2(fftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
% 原书 Listing 2.6：定标正确的二维逆 DFT（乘 (N*delta_f)^2）
N = size(G, 1);
g = fftshift(ifft2(fftshift(G)))*(N*delta_f)^2;
end

function out = rollflip(U)
% 偶数网格居中坐标反转：(Mf)(r) = f(-r)。numpy: roll(flip(U),(1,1))
out = circshift(rot90(U, 2), [1 1]);
end

function [x2, y2, Uout] = one_step_prop(Uin, wvl, d1, Dz)
% 书 Listing 6.1：单步 Fresnel 传播（式 6.5；省略全局相位 e^{ikDz}，原书约定）
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
d2 = wvl*Dz/(N*d1);                            % 式 (6.16)
[x2, y2] = meshgrid(n*d2, n*d2);
Uout = 1/(1i*wvl*Dz) ...
    .* exp(1i*kk/(2*Dz)*(x2.^2 + y2.^2)) ...
    .* ft2(Uin.*exp(1i*kk/(2*Dz)*(x1.^2 + y1.^2)), d1);
end

function [x2, y2, Uout] = two_step_prop(Uin, wvl, d1, d2, Dz)
% 书 Listing 6.3："-"支（Dz1 = Dz/(1-m)）。注意勘误 b：负距离步的中间场
% 坐标对应翻转，m != 1 时输出为真值镜像（本函数忠实复现书印行为）。
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
m = d2/d1;
Dz1 = Dz/(1 - m);
d1a = wvl*abs(Dz1)/(N*d1);                     % 式 (6.19)
[x1a, y1a] = meshgrid(n*d1a, n*d1a);
Uitm = 1/(1i*wvl*Dz1) ...
    .* exp(1i*kk/(2*Dz1)*(x1a.^2 + y1a.^2)) ...
    .* ft2(Uin.*exp(1i*kk/(2*Dz1)*(x1.^2 + y1.^2)), d1);
Dz2 = Dz - Dz1;
[x2, y2] = meshgrid(n*d2, n*d2);
Uout = 1/(1i*wvl*Dz2) ...
    .* exp(1i*kk/(2*Dz2)*(x2.^2 + y2.^2)) ...
    .* ft2(Uitm.*exp(1i*kk/(2*Dz2)*(x1a.^2 + y1a.^2)), d1a);
end

function [x2, y2, Uout] = two_step_prop_fixed(Uin, wvl, d1, d2, Dz)
% "-"支修正版：负距离步的输出按居中反转重排（消除镜像）。
% 数值事实（Python 实测）：修正后与 ang_spec 逐位一致（1e-14，全部测试的 m）。
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
m = d2/d1;
Dz1 = Dz/(1 - m);
d1a = wvl*abs(Dz1)/(N*d1);
[x1a, y1a] = meshgrid(n*d1a, n*d1a);
Uitm = 1/(1i*wvl*Dz1) ...
    .* exp(1i*kk/(2*Dz1)*(x1a.^2 + y1a.^2)) ...
    .* ft2(Uin.*exp(1i*kk/(2*Dz1)*(x1.^2 + y1.^2)), d1);
if Dz1 < 0
    Uitm = rollflip(Uitm);                     % 修正：负距离步输出镜像翻转
end
Dz2 = Dz - Dz1;
[x2, y2] = meshgrid(n*d2, n*d2);
Uout = 1/(1i*wvl*Dz2) ...
    .* exp(1i*kk/(2*Dz2)*(x2.^2 + y2.^2)) ...
    .* ft2(Uitm.*exp(1i*kk/(2*Dz2)*(x1a.^2 + y1a.^2)), d1a);
if Dz2 < 0
    Uout = rollflip(Uout);                     % m<1 时第二步为负距离
end
end

function [x2, y2, Uout] = two_step_prop_plus(Uin, wvl, d1, d2, Dz)
% "+"支（图 6.3 几何）：Dz1 = Dz/(1+m) > 0、Dz2 = m*Dz/(1+m) > 0（m>0）。
N = size(Uin, 1);
kk = 2*pi/wvl;
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
m = d2/d1;
Dz1 = Dz/(1 + m);
d1a = wvl*abs(Dz1)/(N*d1);
[x1a, y1a] = meshgrid(n*d1a, n*d1a);
Uitm = 1/(1i*wvl*Dz1) ...
    .* exp(1i*kk/(2*Dz1)*(x1a.^2 + y1a.^2)) ...
    .* ft2(Uin.*exp(1i*kk/(2*Dz1)*(x1.^2 + y1.^2)), d1);
Dz2 = Dz - Dz1;
[x2, y2] = meshgrid(n*d2, n*d2);
Uout = 1/(1i*wvl*Dz2) ...
    .* exp(1i*kk/(2*Dz2)*(x2.^2 + y2.^2)) ...
    .* ft2(Uitm.*exp(1i*kk/(2*Dz2)*(x1a.^2 + y1a.^2)), d1a);
end

function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
% 书 Listing 6.5：角谱传播（式 6.50 的 m 分支）：
%   Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin / m, d1), df1)
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

function [x2, y2, Uout] = ang_spec_mbranch(Uin, wvl, d1, d2, Dz, corrected)
% 式 (6.66)（m' = -m 分支）的算子链实现。
% corrected=true（= 式 (6.67) 顶符号，正确代换）：
%   U = Q[(m+1)/(m*Dz)] F^-1[f1,-r2/m] Q2[+Dz/m] F Q[(1+m)/Dz] (-1/m)
% corrected=false（书印 6.66）：Q 系数 (1-m)/(m*Dz)、-(1-m)/Dz，F^-1[f1,+r2/m]。
% 实输入场下书印链输出 = -conj(真值)（逐点解析恒等，Python 实测 4e-16）。
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
if corrected
    Q1 = exp(1i*kk/2*(1 + m)/Dz*r1sq);
    Q3 = exp(1i*kk/2*(m + 1)/(m*Dz)*r2sq);
else
    Q1 = exp(-1i*kk/2*(1 - m)/Dz*r1sq);        % 书印 Q[-(1-m)/Dz]
    Q3 = exp(1i*kk/2*(1 - m)/(m*Dz)*r2sq);     % 书印 Q[-(m-1)/(m*Dz)]
end
Q2 = exp(1i*pi^2*2*Dz/m/kk*fsq);               % Q2[+Dz/m]
W = (-1/m)*ift2(Q2 .* ft2(Q1 .* Uin, d1), df1);
if corrected
    W = rollflip(W);                           % F^-1[f1, -r2/m]
end
Uout = Q3 .* W;
end

function [x2, y2, Uout] = ang_spec_propABCD(Uin, wvl, d1, d2, ABCD)
% 书 Listing 6.6：角谱 + ABCD（式 6.79-6.81 的自由 m 版本，笔记补出算子链）：
%   Uout = Q[(D-1/m)/B] * ift2( Q2[-B/m] * ft2( Q[(A-m)/B] * Uin / m ) )
% C 不显式出现（AD-BC=1 保证 (6.79) 自洽）。
% 原书签名注释把 wvl 印成 wwl（笔误），此处按实现更正。
N = size(Uin, 1);
n = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
r1sq = x1.^2 + y1.^2;
df1 = 1/(N*d1);
[fX, fY] = meshgrid(n*df1, n*df1);
fsq = fX.^2 + fY.^2;
m = d2/d1;
[x2, y2] = meshgrid(n*d2, n*d2);
r2sq = x2.^2 + y2.^2;
A = ABCD(1, 1);  B = ABCD(1, 2);  D = ABCD(2, 2);
Q1 = exp(1i*pi/(wvl*B)*(A - m)*r1sq);
Q2 = exp(-1i*pi*wvl*B/m*fsq);
Q3 = exp(1i*pi/(wvl*B)*(D - 1/m)*r2sq);
Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin/m, d1), df1);
end

function U = fresnel_prop_square_ap(x2, y2, D, wvl, Dz)
% 边长 D 方孔的 Fresnel 衍射解析解（式 6.1 可分离积分，省略 e^{ikDz}）：
%   U = (1/(i*lambda*Dz)) I(x2) I(y2)，
%   I(x) = sqrt(lambda Dz/2) * {[C(t2)-C(t1)] + i[S(t2)-S(t1)]}，
%   t = sqrt(2/(lambda Dz)) * u
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
% Fresnel 积分 S(t)=int_0^t sin(pi u^2/2)du、C(t)=int_0^t cos(pi u^2/2)du。
% MATLAB 基础版无内置（符号工具箱除外）；这里用 integral() 数值积分
% （Python 版用 scipy.special.fresnel）。t 为数组，逐元素求值。
S = arrayfun(@(tt) integral(@(u) sin(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
C = arrayfun(@(tt) integral(@(u) cos(pi*u.^2/2), 0, tt, ...
    'AbsTol', 1e-10, 'RelTol', 1e-8), t);
end

function out = rect(x)
% rect(x)：|x|<1/2 取 1，|x|=1/2 取 1/2，其余 0
out = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end

function out = circ_(x, y, R)
% 半径 R 圆盘（半径参数；原书 circ(x,y,D) 是直径参数，注意换算）
out = double(hypot(x, y) <= R);
end
