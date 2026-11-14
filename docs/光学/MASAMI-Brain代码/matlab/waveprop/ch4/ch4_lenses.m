%% ===== 书 Schmidt PM199 第 4 章 =====
%  透镜的傅里叶变换性质（原书 4.2 节，Listings 4.3-4.5 与图 4.2-4.4）
%  与 python/waveprop/ch4/ch4_lenses.py 数值行为一致（以 Python 实跑为准）。
%
%  焦距 f_l 的会聚薄镜给场乘相位 exp(-i k/(2 f_l)(x^2+y^2))（式 4.7）。
%  物片放三个位置、传到后焦面（Delta_z = f_l）时，透镜相位与菲涅耳积分
%  （式 4.6）的二次相位相消，剩下的都是"prefactor x 相位因子 x 一个 FT"：
%      物贴透镜      式 (4.10)：1/(i lambda f_l) x exp(i k/(2 f_l) r2^2)
%      物在透镜前 d  式 (4.12)：相位 exp(i k/(2 f_l)(1-d/f_l) r2^2)，
%                    严格形式带移位光瞳 P(x1 + d/f_l x2)（渐晕）——
%                    Listing 4.4 未实现该项（相当于 P=1）
%      物在透镜后 d  式 (4.16)：(f_l/d)/(i lambda d) x exp(i k/(2 d) r2^2)，
%                    FT 标度 lambda*d（几何光学近似，要求 d << f_l）
%  特例：d=0 退回物贴透镜；d=f_l（物在前焦面）相位因子=1，精确 FT。
%
%  原书勘误（本脚本按修正签名实现）：
%      Listing 4.5 第 1-2 行签名 lens_behind_ft(Uin, wvl, d1, f) 少了 d，
%      函数体却用到 d；第 4 行注释 (Uin, wvl, d1, d, f) 才是对的。
%      另外 Listing 4.4/4.5 的返回值顺序 [x2 y2 Uout] 与 4.1/4.3 的
%      [Uout x2 y2] 不一致（本脚本统一为 [Uout x2 y2]）。
%
%  验证（编号与 Python 一致）：
%   (1) d=0 时 lens_in_front_ft == lens_against_ft（逐点相等）
%   (2) d=f_l 退化为精确 FT；与贴透镜差恰好 exp(i k/(2 f_l) r2^2)
%   (3) 独立路线：直接按式 (4.6) 二维求和（不用 FFT）vs FT 路线
%   (4) 能量守恒：sum|U2|^2 d2^2 == 物面积
%   (5) 物在透镜后：主瓣全宽 2*lambda*d/D 与峰值 (f_l/d)^2*(D^2/(lambda d))^2
%   (6) 渐晕：移位光瞳 vs Listing 4.4 —— 主瓣内一致，主瓣外 sinc 零点被填掉
clear; clc; close all;

%% 演示参数：0.5 mm 方孔物片，f_l = 1 m，lambda = 1 um
N = 256;  L = 2.5e-3;  d1 = L/N;  D = 0.5e-3;
wvl = 1e-6;  k = 2*pi/wvl;  fl = 1.0;
[x1, y1] = meshgrid((-N/2:N/2-1)*d1);
tA = rect_fun(x1/D) .* rect_fun(y1/D);
area = D^2;

%% 验证 1：d=0 等价
[U_ag, x2, y2] = lens_against_ft(tA, wvl, d1, fl);
U_bf0 = lens_in_front_ft(tA, wvl, d1, fl, 0.0);
fprintf('验证 1：max|U_ag - U_bf(d=0)| = %.2e\n', max(abs(U_ag(:)-U_bf0(:))));

%% 验证 2：d=f_l 精确 FT
U_bf1 = lens_in_front_ft(tA, wvl, d1, fl, fl);
U_exact = 1/(1i*wvl*fl) * ft2(tA, d1);
fprintf('验证 2：max|angle(U_bf(d=f_l)/U_exact)| = %.2e rad\n', ...
    max(abs(angle(U_bf1 ./ U_exact))));
U_ag2 = U_exact .* exp(1i*k/(2*fl)*(x2.^2+y2.^2));
fprintf('        max|U_ag - exp(ik/2f_l r2^2)*U_exact| = %.2e\n', ...
    max(abs(U_ag(:)-U_ag2(:))));
fprintf('        贴透镜时 r2=1 cm 处焦面相位 = %.1f rad\n', k/(2*fl)*0.01^2);

%% 验证 3：两条路线（直接积分 vs FT）
U1 = tA .* exp(-1i*k/(2*fl)*(x1.^2+y1.^2));   % 片后场 = tA x 透镜相位
xs = x2(N/2+1, :);                            % y2=0 行的 x2
U_dir = fresnel_direct_slice(U1, wvl, d1, fl, xs, 0.0);
rel = max(abs(abs(U_dir) - abs(U_ag(N/2+1,:)))) / max(abs(U_ag(N/2+1,:)));
fprintf('验证 3：直接积分 vs FT：max 相对差 = %.2e\n', rel);

%% 验证 4：能量守恒
d2 = wvl*fl/L;
P_num = sum(abs(U_ag(:)).^2) * d2^2;
fprintf('验证 4：sum|U2|^2 d2^2 = %.5e vs 物面积 = %.5e（相对偏差 %.2e）\n', ...
    P_num, area, abs(P_num-area)/area);
fprintf('        轴上辐照度 = %.6f W/m^2，理论 (D^2/(lambda f_l))^2 = %.6f\n', ...
    abs(U_ag(N/2+1,N/2+1))^2, (area/(wvl*fl))^2);

%% 验证 5：物在透镜后——宽度与峰值随 d
for d = [0.5*fl, 0.25*fl]
    U_bh = lens_behind_ft(tA, wvl, d1, d, fl);
    row = abs(U_bh(N/2+1, :)).^2;
    [pk, i_pk] = max(row);
    x2b = wvl*d*((-N/2:N/2-1)/(N*d1));        % 该函数的观察面坐标
    right = i_pk + find(row(i_pk:end) < 0.01*pk, 1) - 1;
    mainlobe = 2*(x2b(right) - x2b(i_pk));
    peak_th = (fl/d)^2 * (area/(wvl*d))^2;
    fprintf('验证 5：d=%.2f m：峰值 %.4e vs 理论 %.4e；主瓣全宽 %.3e vs 理论 %.3e m\n', ...
        d, pk, peak_th, mainlobe, 2*wvl*d/D);
end

%% 验证 6：渐晕——式 (4.11) 的移位光瞳 vs Listing 4.4
d = 0.5*fl;  Dp = 2e-3;
P0 = double(sqrt(x1.^2+y1.^2) <= Dp/2);
U_no = lens_in_front_ft(tA .* P0, wvl, d1, fl, d);
inside = tA ~= 0;
xa = x1(inside);  ya = y1(inside);
m6 = abs(xs) <= 2.5e-3;
xs6 = xs(m6);
U_strict = zeros(size(xs6));
for j = 1:numel(xs6)
    sh = d/fl*xs6(j);
    Pv = double(sqrt((xa+sh).^2+ya.^2) <= Dp/2);   % 移位光瞳
    U_strict(j) = sum(tA(inside) .* Pv ...
        .* exp(-2i*pi/(wvl*fl)*xs6(j)*xa)) * d1^2/(1i*wvl*fl) ...
        * exp(1i*k/(2*fl)*(1-d/fl)*xs6(j)^2);
end
U_no_sl = abs(U_no(N/2+1, m6));
for xv = [0, 1e-3, 1.5e-3, 2e-3]
    [~, j] = min(abs(xs6 - xv));
    fprintf('验证 6：x2 = %.1f mm：|U_严格|/|U_Listing4.4| = %.3f\n', ...
        xv*1e3, abs(U_strict(j))/U_no_sl(j));
end
fprintf('        （主瓣内一致；主瓣外移位光瞳填掉 sinc 零点。'
    'Listing 4.4 只在 P=1 或光瞳随物片移动时严格成立）\n');

%% 绘图：图 4-7
U_bf05 = lens_in_front_ft(tA, wvl, d1, fl, 0.5*fl);
U_bh = lens_behind_ft(tA, wvl, d1, 0.5*fl, fl);
m5 = abs(x2(N/2+1,:)) <= 2.5e-3;
figure('Position',[40 40 1150 420]);
subplot(1,2,1);
plot(x2(N/2+1,m5)*1e3, abs(U_ag(N/2+1,m5)).^2, 'k-', 'LineWidth', 1.5); hold on;
plot(x2(N/2+1,m5)*1e3, abs(U_bf05(N/2+1,m5)).^2, 'b.', 'MarkerSize', 6);
m5b = abs(x2(N/2+1,:)) <= 2.5e-3;   % 物后的坐标轴标度是 lambda*d
x2b5 = wvl*0.5*((-N/2:N/2-1)/(N*d1));
plot(x2b5(m5b)*1e3, abs(U_bh(N/2+1,m5b)).^2, 'r--', 'LineWidth', 1.3);
xlabel('x_2 [mm]'); ylabel('辐照度 [W/m^2]'); grid on;
title('(a) 焦面辐照度：物前位置不影响图样，物后改变标度');
legend({'物贴透镜（=物前 d=0）', '物前 d=0.5f_l', '物后 d=0.5f_l（标度 \lambda d）'});
subplot(1,2,2);
m7 = abs(x2(N/2+1,:)) <= 2e-3;
x7 = x2(N/2+1, m7);
xfine = linspace(-2e-3, 2e-3, 2001);
cols = {'k', [0.2 0.4 0.85], [0 0.6 0]};
labs = {'物贴透镜（曲率 1/f_l）', '物前 d=0.5f_l（0.5/f_l）', ...
        '物前 d=f_l（精确 FT，零相位）'};
ds = [0, 0.5*fl, fl];
for q = 1:3
    phi = angle(exp(1i*k/(2*fl)*(1-ds(q)/fl)*xfine.^2));
    plot(xfine*1e3, phi, '-', 'Color', cols{q}, 'LineWidth', 1.1); hold on;
end
plot(x7*1e3, angle(U_ag(N/2+1,m7)), 'k.', 'MarkerSize', 7);
plot(x7*1e3, angle(U_bf05(N/2+1,m7)), 'b.', 'MarkerSize', 7);
xlabel('x_2 [mm]'); ylabel('相位 [rad]'); grid on;
title('(b) 焦面相位记得物的位置（点=数值，线=理论包裹相位）');
legend(labs, 'FontSize', 8, 'Location', 'south');

%% ===== 本地函数 =====
function Uout = lens_against_ft(Uin, wvl, d1, f)
% 原书 Listing 4.3：物贴透镜 -> 后焦面（式 4.10）
N = size(Uin, 1);
k = 2*pi/wvl;
fX = (-N/2 : N/2-1) / (N*d1);
[x2, y2] = meshgrid(wvl * f * fX);        % 观察面坐标标度 lambda*f_l
Uout = exp(1i*k/(2*f)*(x2.^2+y2.^2)) / (1i*wvl*f) .* ft2(Uin, d1);
end

function Uout = lens_in_front_ft(Uin, wvl, d1, f, d)
% 原书 Listing 4.4：物在透镜前 d -> 后焦面（式 4.12 取 P=1）。
% 原书返回值顺序 [x2 y2 Uout]，这里统一为 Uout。
N = size(Uin, 1);
k = 2*pi/wvl;
fX = (-N/2 : N/2-1) / (N*d1);
[x2, y2] = meshgrid(wvl * f * fX);
Uout = 1/(1i*wvl*f) ...
    .* exp(1i*k/(2*f)*(1-d/f)*(x2.^2+y2.^2)) ...
    .* ft2(Uin, d1);
end

function Uout = lens_behind_ft(Uin, wvl, d1, d, f)
% 原书 Listing 4.5：物在透镜后 -> 焦面（式 4.16）。
% 勘误：原书签名缺 d，本签名 (Uin, wvl, d1, d, f) 按其函数体补正。
N = size(Uin, 1);
k = 2*pi/wvl;
fX = (-N/2 : N/2-1) / (N*d1);
[x2, y2] = meshgrid(wvl * d * fX);        % 标度是 lambda*d ！
Uout = f/d * 1/(1i*wvl*d) ...
    .* exp(1i*k/(2*d)*(x2.^2+y2.^2)) ...
    .* ft2(Uin, d1);
end

function G = ft2(g, delta)
% 原书 Listing 2.5
G = fftshift(fft2(fftshift(g))) * delta^2;
end

function s = rect_fun(x)
% rect(x)：|x|<1/2 取 1，|x|=1/2 取 1/2，其余 0
s = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end

function U = fresnel_direct_slice(U1, wvl, d1, Dz, x2_vals, y2)
% 独立路线：按菲涅耳衍射积分（式 4.6）直接二维求和（不用 FFT），
% 只算 y=y2 一行。分块避免 N^2 x N^2 的大矩阵。
N = size(U1, 1);
k = 2*pi/wvl;
xa = (-N/2:N/2-1)*d1;
[X1, Y1] = meshgrid(xa, xa);
src = (U1 .* exp(1i*k/(2*Dz)*(X1.^2+Y1.^2)));   % 源面二次相位
src = reshape(src, 1, []);
Xv = reshape(X1, 1, []);  Yv = reshape(Y1, 1, []);
U = zeros(size(x2_vals));
blk = 32;
for i0 = 1:blk:numel(x2_vals)
    idx = i0:min(i0+blk-1, numel(x2_vals));
    K = exp(-2i*pi/(wvl*Dz) * (x2_vals(idx).' * Xv + y2*Yv));
    U(idx) = K * src.';
end
U = U * d1^2 / (1i*wvl*Dz);
end
