%% ===== 书 Schmidt PM199 第 4 章 =====
%  原书 4.3 节习题 1-5：夫琅禾费衍射与透镜焦面的五个数值/解析对比
%  与 python/waveprop/ch4/ch4_problems.py 数值行为一致（以 Python 实跑为准）。
%
%  习题（网格默认按 4.1 节例题：N=512、L=7.5 mm、Dz=20 m、lambda=1 um）：
%   1) 1 mm 方孔：U = D^2 sinc(D fx) sinc(D fy)，fx = x2/(lambda*Dz)；
%   2) 双缝：两个 1 mm 方孔、边缘间隙 0.5 mm（中心距 1.5 mm）：
%      2 D^2 sinc(D fx) sinc(D fy) cos(2 pi fx x0)，x0 = 0.75 mm；
%   3) 振幅光栅 f0=10/D：D^2 sinc(D fy)[1/2 sinc(D fx)
%      + 1/4 sinc(D(fx-f0)) + 1/4 sinc(D(fx+f0))]；
%   4) 相位光栅 f0=10/D：e^{i 2 pi cos(2 pi f0 x)} 的 Jacobi-Anger 展开
%      给出衍射级，第 n 级峰值强度 J_n^2(2 pi)(D^2/(lambda Dz))^2，
%      sum J_n^2(2pi)=1。带宽账本：J_4..J_6 都不小（0.32/0.37/0.28），
%      N=512 混叠 69%、N=1024 混叠 6.2% -> 本题取 N=2048；
%   5) 高斯束聚焦（w=2 cm 束腰在透镜、f_l=1 m、无限孔径）：
%      解析焦面 = pi w^2/(i lambda f_l) exp(-(pi w r2/(lambda f_l))^2)，
%      强度 1/e^2 半径 = lambda f_l/(pi w) = 15.9 um。
%      网格单独设计：L=0.4 m、N=1024（束要大窗、焦斑要细 d2）。
%
%  MATLAB 注意：besselj(n,2*pi) 对应 scipy jv(n,2*pi)；
%  sinc 用归一化定义（本地函数，避免依赖工具箱）。需要 R2018b+。
clear; clc; close all;

%% 公共网格（4.1 节例题）
N = 512;  L = 7.5e-3;  d1 = L/N;  wvl = 1e-6;  k = 2*pi/wvl;  Dz = 20;
[x1, y1] = meshgrid((-N/2:N/2-1)*d1);
D = 1e-3;  f0 = 10/D;

%% 习题 1：方孔
Uin = rect_fun(x1/D) .* rect_fun(y1/D);
[Uout, x2, y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
A1 = D^2 * sinc_fun(D*x2/(wvl*Dz)) .* sinc_fun(D*y2/(wvl*Dz));
I_th = abs(exp(1i*k/(2*Dz)*(x2.^2+y2.^2))/(1i*wvl*Dz) .* A1).^2;
I_num = abs(Uout).^2;
report('习题 1', I_num, I_th);
res1 = {x2(N/2+1,:), I_num(N/2+1,:), I_th(N/2+1,:)};

%% 习题 2：双缝（边缘间隙 0.5 mm）
x0 = 0.75e-3;
Uin = (rect_fun((x1-x0)/D) + rect_fun((x1+x0)/D)) .* rect_fun(y1/D);
[Uout, x2, y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
A2 = 2*D^2 * sinc_fun(D*x2/(wvl*Dz)) .* sinc_fun(D*y2/(wvl*Dz)) ...
    .* cos(2*pi*x2/(wvl*Dz)*x0);
I_th = abs(exp(1i*k/(2*Dz)*(x2.^2+y2.^2))/(1i*wvl*Dz) .* A2).^2;
I_num = abs(Uout).^2;
report('习题 2', I_num, I_th);
res2 = {x2(N/2+1,:), I_num(N/2+1,:), I_th(N/2+1,:)};

%% 习题 3：振幅光栅
Uin = 0.5*(1+cos(2*pi*f0*x1)) .* rect_fun(x1/D) .* rect_fun(y1/D);
[Uout, x2, y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
fxg = x2/(wvl*Dz);
A3 = D^2 * sinc_fun(D*y2/(wvl*Dz)) .* ...
    (0.5*sinc_fun(D*fxg) + 0.25*sinc_fun(D*(fxg-f0)) + 0.25*sinc_fun(D*(fxg+f0)));
I_th = abs(exp(1i*k/(2*Dz)*(x2.^2+y2.^2))/(1i*wvl*Dz) .* A3).^2;
I_num = abs(Uout).^2;
report('习题 3', I_num, I_th);
res3 = {x2(N/2+1,:), I_num(N/2+1,:), I_th(N/2+1,:)};

%% 习题 4：相位光栅（N=2048）
N4 = 2048;  d14 = L/N4;
[x1b, y1b] = meshgrid((-N4/2:N4/2-1)*d14);
Uin = exp(2i*pi*cos(2*pi*f0*x1b)) .* rect_fun(x1b/D) .* rect_fun(y1b/D);
[Uout, x2b, y2b] = fraunhofer_prop(Uin, wvl, d14, Dz);
I_num = abs(Uout).^2;
A4 = zeros(size(x2b));
for n_ord = -8:8
    A4 = A4 + 1i^n_ord * besselj(n_ord, 2*pi) * D^2 ...
        .* sinc_fun(D*(x2b/(wvl*Dz)-n_ord*f0)) .* sinc_fun(D*y2b/(wvl*Dz));
end
I_th = abs(exp(1i*k/(2*Dz)*(x2b.^2+y2b.^2))/(1i*wvl*Dz) .* A4).^2;
report('习题 4', I_num, I_th);
scale_pk = (D^2/(wvl*Dz))^2;
row4 = I_num(N4/2+1, :);
fx_row = x2b(N4/2+1, :) / (wvl*Dz);
for n_ord = 0:3
    [~, j] = min(abs(fx_row - n_ord*f0));
    fprintf('    n=%d：数值 %.4e vs 理论 %.4e（J_%d(2pi)=%.4f）\n', ...
        n_ord, row4(j), besselj(n_ord,2*pi)^2*scale_pk, n_ord, ...
        besselj(n_ord,2*pi));
end
s_all = 0;
for n_ord = -60:60
    s_all = s_all + besselj(abs(n_ord), 2*pi)^2;
end
for q = [512, 1024, 2048]
    nb = floor(1/(2*(L/q)*f0));
    s_in = 0;
    for n_ord = -nb:nb
        s_in = s_in + besselj(abs(n_ord), 2*pi)^2;
    end
    fprintf('    N=%d：频带 ±%.0f m^-1 装到 |n|<=%d，能量占比 %.1f%%（混叠 %.2f%%）\n', ...
        q, 1/(2*(L/q)), nb, 100*s_in/s_all, 100*(1-s_in/s_all));
end
P4 = sum(I_num(:)) * (wvl*Dz/L)^2;
fprintf('    能量：sum|U2|^2 d2^2 = %.4e vs D^2 = %.4e（相对偏差 %.2e）\n', ...
    P4, D^2, abs(P4-D^2)/D^2);
res4 = {x2b(N4/2+1,:), row4, I_th(N4/2+1,:)};

%% 习题 5：高斯束聚焦
w = 2e-2;  fl = 1.0;  L5 = 0.4;  N5 = 1024;  d15 = L5/N5;
[x1c, y1c] = meshgrid((-N5/2:N5/2-1)*d15);
Uin5 = exp(-(x1c.^2+y1c.^2)/w^2);
[Uout5, x2c, y2c] = lens_against_ft(Uin5, wvl, d15, fl);
r2 = sqrt(x2c.^2+y2c.^2);
U_th5 = pi*w^2/(1i*wvl*fl) .* exp(-(pi*w*r2/(wvl*fl)).^2);
I5n = abs(Uout5).^2;  I5a = abs(U_th5).^2;
shape5 = abs(I5n/max(I5n(:)) - I5a/max(I5a(:)));
fprintf('习题 5：采样 d1=%.3e m（每 w 跨 %.0f 样本），d2=%.1e m\n', ...
    d15, w/d15, wvl*fl/L5);
fprintf('        峰值 %.4e vs 理论 %.4e；形状误差 max=%.2e\n', ...
    max(I5n(:)), (pi*w^2/(wvl*fl))^2, max(shape5(:)));
row5 = I5n(N5/2+1, :);  x5row = x2c(N5/2+1, :);
half = max(row5)/exp(2);
[~, i_pk] = max(row5);
right = i_pk + find(row5(i_pk:end) < half, 1) - 1;
xcross = x5row(right) - (row5(right)-half)/(row5(right)-row5(right-1)) ...
    * (x5row(right)-x5row(right-1));
fprintf('        1/e^2 半径：数值 %.4e m vs 理论 lambda f_l/(pi w) = %.4e m\n', ...
    xcross - x5row(i_pk), wvl*fl/(pi*w));
P5 = sum(I5n(:)) * (wvl*fl/L5)^2;
fprintf('        能量：%.4e vs pi w^2/2 = %.4e（相对偏差 %.2e）\n', ...
    P5, pi*w^2/2, abs(P5-pi*w^2/2)/P5);

%% 绘图：图 4-8 与图 4-9
res_all = {res1, res2, res3, res4};
titles = {'习题 1：1 mm 方孔', '习题 2：双缝（间隙 0.5 mm）', ...
    '习题 3：振幅光栅 f_0=10/D', '习题 4：相位光栅 f_0=10/D（N=2048）'};
halves = [0.1, 0.1, 0.3, 0.7];
figure('Position',[30 30 1150 760]);
for q = 1:4
    r = res_all{q};
    m = abs(r{1}) <= halves(q);
    subplot(2,2,q);
    plot(r{1}(m)*1e3, r{3}(m), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); hold on;
    plot(r{1}(m)*1e3, r{2}(m), 'b.', 'MarkerSize', 3);
    xlabel('x_2 [mm]'); ylabel('辐照度 [W/m^2]'); grid on;
    title(titles{q}, 'FontSize', 11); legend({'解析', '数值'});
end

figure('Position',[40 40 1150 420]);
subplot(1,2,1);
msrc = abs(x1c(1,:)) <= 0.05;
plot(x1c(1,msrc), Uin5(N5/2+1,msrc), 'k-', 'LineWidth', 1.3); hold on;
xline(-w, 'r:'); xline(w, 'r:');
xlabel('x_1 [m]'); ylabel('振幅'); grid on;
title('(a) 透镜面振幅（红线 ±w，窗 L=0.4 m 截断 e^{-100}）');
subplot(1,2,2);
mf = abs(x5row) <= 6e-5;
plot(x5row(mf)*1e6, I5a(mf), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6); hold on;
plot(x5row(mf)*1e6, row5(mf), 'b.', 'MarkerSize', 5);
xlabel('x_2 [\mum]'); ylabel('辐照度 [W/m^2]'); grid on;
title('(b) 焦面光斑（1/e^2 半径 15.9 \mum）'); legend({'解析', '数值'});

%% ===== 本地函数 =====
function [Uout, x2, y2] = fraunhofer_prop(Uin, wvl, d1, Dz)
% 原书 Listing 4.1
N = size(Uin, 1);
k = 2*pi/wvl;
fX = (-N/2 : N/2-1) / (N*d1);
[x2, y2] = meshgrid(wvl * Dz * fX);
Uout = exp(1i*k/(2*Dz)*(x2.^2+y2.^2)) / (1i*wvl*Dz) .* ft2(Uin, d1);
end

function Uout = lens_against_ft(Uin, wvl, d1, f)
% 原书 Listing 4.3
N = size(Uin, 1);
k = 2*pi/wvl;
fX = (-N/2 : N/2-1) / (N*d1);
[x2, y2] = meshgrid(wvl * f * fX);
Uout = exp(1i*k/(2*f)*(x2.^2+y2.^2)) / (1i*wvl*f) .* ft2(Uin, d1);
end

function G = ft2(g, delta)
% 原书 Listing 2.5
G = fftshift(fft2(fftshift(g))) * delta^2;
end

function s = rect_fun(x)
s = double(abs(x) < 1/2) + 0.5*double(abs(x) == 1/2);
end

function s = sinc_fun(x)
% 归一化 sinc：sin(pi x)/(pi x)（避免依赖 Signal Processing Toolbox）
s = ones(size(x));
nz = x ~= 0;
s(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end

function report(tag, I_num, I_th)
% 标度差 + 归一化形状误差（与 Python 版 err_report 一致）
shape = abs(I_num/max(I_num(:)) - I_th/max(I_th(:)));
fprintf('%s：标度差 %+.2f%%，归一化形状误差 max %.2e\n', ...
    tag, 100*(max(I_num(:))/max(I_th(:))-1), max(shape(:)));
end
