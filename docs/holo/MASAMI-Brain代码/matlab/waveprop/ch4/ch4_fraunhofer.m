%% ===== 书 Schmidt PM199 第 4 章 =====
%  夫琅禾费衍射（原书 4.1 节，Listings 4.1-4.2 与图 4.1）
%  与 python/waveprop/ch4/ch4_fraunhofer.py 数值行为一致（本机无 MATLAB 时
%  以 Python 实跑结果为准，见笔记页面打印输出）。
%
%  式 (4.3)：夫琅禾费衍射 = 一次定标正确的 ft2，观察面坐标按
%  x2 = lambda*Dz*fx1 线性映射（式 4.4）。远场判据（式 4.2）：
%  Delta_z > 2*D^2/lambda（本例 = 2000 m，例题只用 20 m —— 该判据保守，
%  真正起作用的是孔径上二次相位 < ~pi，即 Delta_z >> D^2/(4*lambda) = 0.25 m）。
%
%  验证（编号与 Python 一致）：
%   (0) jinc 约定：jinc(u) = 2*J1(pi*u)/(pi*u)，jinc(0)=1（第一零点 u=1.22）
%   (1) 书例复现：轴上辐照度、传感器区域误差（分"标度误差/形状误差"）
%   (2) 能量守恒：sum|U2|^2*d2^2 == 孔径面积（离散 Parseval，混叠也守恒）
%   (3) 采样不足：N=64 -> 观察窗变窄，图样卷回窗内
%   (4) 两条路线：直接按菲涅耳积分 (4.6) 二维求和 vs FT 路线
%
%  MATLAB 注意：数组 1 起，居中排列的零坐标样本在 N/2+1 号；
%  meshgrid 与 np.meshgrid 默认 'xy' 行为一致。需要 R2018b+（本地函数与
%  xline；旧版把 xline 换成 plot 两点即可）。
clear; clc; close all;

%% 参数（原书 Listing 4.2 的 example_fraunhofer_circ.m）
N = 512;            % number of grid points per side
L = 7.5e-3;         % total size of the grid [m]
d1 = L / N;         % source-plane grid spacing [m]
D = 1e-3;           % diameter of the aperture [m]
wvl = 1e-6;         % optical wavelength [m]
k = 2*pi / wvl;
Dz = 20;            % propagation distance [m]
area = pi*D^2/4;    % 圆孔面积 = 谱的直流值

[x1, y1] = meshgrid((-N/2 : N/2-1) * d1);
Uin = circ_fun(x1, y1, D);

%% 验证 0：jinc 约定核对
rho = linspace(1e-6, 400, 400001);
diff0 = max(abs((D/2)*besselj(1, pi*D*rho)./rho - area .* jinc_s(D*rho)));
fprintf('验证 0：max|base_jinc(rho) - 面积*jinc_s(D*rho)| = %.2e\n', diff0);
fprintf('        jinc_s 第一零点 u=%.4f -> 爱里环半径 %.4f m\n', ...
    3.8317/pi, 3.8317/pi*wvl*Dz/D);

%% 验证 1 + 2：书例复现与能量守恒
[Uout, x2, y2] = fraunhofer_prop(Uin, wvl, d1, Dz);
d2 = wvl*Dz/L;  x2max = wvl*Dz/(2*d1);
fprintf('验证 1：d1=%.3e m, d2=%.4f m, 观察窗 ±%.3f m\n', d1, d2, x2max);
U_th = exp(1i*k/(2*Dz)*(x2.^2+y2.^2)) / (1i*wvl*Dz) ...
    * area * jinc_s(D*sqrt(x2.^2+y2.^2)/(wvl*Dz));
I_num = abs(Uout).^2;  I_ana = abs(U_th).^2;
peak = (area/(wvl*Dz))^2;
fprintf('        轴上辐照度：数值 %.6f vs 解析 %.6f mW/m^2（高 %.2f%%）\n', ...
    I_num(N/2+1,N/2+1)*1e3, peak*1e3, 100*(I_num(N/2+1,N/2+1)/peak-1));
shape = abs(I_num/max(I_num(:)) - I_ana/max(I_ana(:)));
sensor = abs(x2) <= 0.2;
[~, imax] = max(shape(N/2,:));
fprintf('        传感器区域：标度误差 %.2e，归一化形状误差 %.2e\n', ...
    max(abs(I_num(:)-I_ana(:)).*(sensor(:)))/peak, max(shape(sensor)));
fprintf('        形状误差最大值 %.2e 在 x2=%.3f m（爱里环陡坡）；', ...
    max(shape(:)), x2(N/2, imax));
fprintf('边缘（|x2|>0.6 m）仅 %.1e\n', max(shape(abs(x2)>0.6)));
P_num = sum(I_num(:)) * d2^2;
fprintf('验证 2：sum|U2|^2*d2^2 = %.6e vs 孔径面积 %.6e（相对偏差 %.2e）\n', ...
    P_num, area, abs(P_num-area)/area);
nf = 1501; xf = linspace(-x2max, x2max, nf);
[Xf, Yf] = meshgrid(xf);
I_fine = (area/(wvl*Dz))^2 * jinc_s(D*sqrt(Xf.^2+Yf.^2)/(wvl*Dz)).^2;
P_ana_win = trapz(xf, trapz(xf, I_fine));
fprintf('        解析场窗内能量占 %.2f%% -> 混叠卷回约 %.2f%%\n', ...
    100*P_ana_win/area, 100*(1-P_ana_win/area));

%% 验证 3：采样不足的混叠（N=64）
N2 = 64; d1b = L/N2;
[x1b, y1b] = meshgrid((-N2/2:N2/2-1)*d1b);
[Uout_b, x2b, y2b] = fraunhofer_prop(circ_fun(x1b, y1b, D), wvl, d1b, Dz);
U_thb = exp(1i*k/(2*Dz)*(x2b.^2+y2b.^2)) / (1i*wvl*Dz) ...
    * area * jinc_s(D*sqrt(x2b.^2+y2b.^2)/(wvl*Dz));
I_numb = abs(Uout_b).^2;  I_anab = abs(U_thb).^2;
shapeb = abs(I_numb/max(I_numb(:)) - I_anab/max(I_anab(:)));
fprintf('验证 3：N=64 观察窗 ±%.4f m，形状误差 max=%.3f（外环错位），', ...
    wvl*Dz/(2*d1b), max(shapeb(:)));
fprintf('标度误差 %.2f%%\n', 100*(max(I_numb(:))/max(I_anab(:))-1));
P_numb = sum(I_numb(:)) * (wvl*Dz/L)^2;
xf_b = linspace(-wvl*Dz/(2*d1b), wvl*Dz/(2*d1b), 2001);
[Xfb, Yfb] = meshgrid(xf_b);
Ifb = (area/(wvl*Dz))^2 * jinc_s(D*sqrt(Xfb.^2+Yfb.^2)/(wvl*Dz)).^2;
P_anab = trapz(xf_b, trapz(xf_b, Ifb));
fprintf('        能量账本：数值仍=%.4e（=面积）；解析窗内=%.4e（%.2f%%）\n', ...
    P_numb, P_anab, 100*P_anab/area);

%% 验证 4：两条路线对照（直接积分 vs FT）
inside = Uin ~= 0;
xa = x1(inside);  ya = y1(inside);
quad = exp(1i*k/(2*Dz)*(xa.^2+ya.^2));    % 被远场近似扔掉的二次相位
xs = x2(N/2, :);                          % y2=0 行
K = exp(-1i*k/Dz * (xs.' * xa));          % (N x Ninside)
U_dir = (K*quad) * d1^2 / (1i*wvl*Dz) .* exp(1i*k/(2*Dz)*xs.^2).';
rel = max(abs(abs(U_dir) - abs(Uout(N/2,:)))) / max(abs(Uout(N/2,:)));
fprintf('验证 4：直接积分 vs FT：max 相对差 = %.2e（孔径边缘二次相位仅 %.3f rad）\n', ...
    rel, k/(2*Dz)*(D/2)^2);

%% 绘图：图 4-5（书例复现）与图 4-6（混叠）
half = round(0.1/d2) + 1;
sl = N/2+1 + (-half:half-1);
figure('Position',[40 40 1350 400]);
subplot(1,3,1);
imagesc(x2(N/2,sl), x2(N/2,sl), I_num(sl, sl)); axis image; colorbar;
caxis([0 max(I_num(:))]); colormap inferno;
xlabel('x_2 [m]'); ylabel('y_2 [m]');
title('(a) 观察面辐照度（±0.1 m）');
subplot(1,3,2);
plot(x2(N/2,sl), I_ana(N/2,sl), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6); hold on;
plot(x2(N/2,sl), I_num(N/2,sl), 'b.', 'MarkerSize', 5);
xlabel('x_2 [m]'); ylabel('辐照度 [W/m^2]'); grid on;
title('(b) y_2=0 切片：数值 vs 解析'); legend({'解析', '数值'});
subplot(1,3,3);
plot(x2(N/2,sl), shape(N/2,sl), 'k.', 'MarkerSize', 4);
xlabel('x_2 [m]'); ylabel('归一化形状误差'); grid on;
title('(c) 残差在爱里环陡坡处最大');

hb = N2/2 + 1 + (-N2/2:N2/2-1);
figure('Position',[40 60 1050 400]);
subplot(1,2,1);
plot(x2b(hb,1), I_anab(hb,:), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6); hold on;
plot(x2b(hb,1), I_numb(hb,:), 'b.', 'MarkerSize', 5);
xline(wvl*Dz/(2*d1b), 'r:'); xline(-wvl*Dz/(2*d1b), 'r:');
xlabel('x_2 [m]'); ylabel('辐照度 [W/m^2]'); grid on;
title('(a) N=64：观察窗变窄，图样卷回窗内'); legend({'解析', '数值（N=64）'});
subplot(1,2,2);
semilogy(abs(x2(N/2,sl)), shape(N/2,sl)+1e-16, 'k.', 'MarkerSize', 3); hold on;
semilogy(abs(x2b(hb,:)), shapeb(hb,:)+1e-16, 'r.', 'MarkerSize', 4);
xlabel('|x_2| [m]'); ylabel('归一化形状误差'); grid on;
title('(b) 同一位置的误差：N=64（红）比 N=512（黑）大两个量级');
legend({'N=512（窗 ±0.68 m）', 'N=64（窗 ±0.085 m）'});

%% ===== 本地函数 =====
function [Uout, x2, y2] = fraunhofer_prop(Uin, wvl, d1, Dz)
% 原书 Listing 4.1：夫琅禾费传播（e^{ikDz} 轴上相位不计算）
N = size(Uin, 1);                 % assume square grid
k = 2*pi / wvl;
fX = (-N/2 : N/2-1) / (N*d1);     % 源面频率轴（升序、居中排列）
[x2, y2] = meshgrid(wvl * Dz * fX);   % 式 (4.4)：x2 = lambda*Dz*fx1
Uout = exp(1i*k/(2*Dz)*(x2.^2+y2.^2)) / (1i*wvl*Dz) .* ft2(Uin, d1);
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）
G = fftshift(fft2(fftshift(g))) * delta^2;
end

function U = circ_fun(X, Y, D)
% 直径 D 的圆孔（原书 circ(x,y,D) 约定：|r| <= D/2 取 1）
U = double(sqrt(X.^2+Y.^2) <= D/2);
end

function out = jinc_s(u)
% 原书解析式用的 jinc：jinc(u) = 2*J1(pi*u)/(pi*u)，jinc(0)=1
u(u == 0) = eps;
out = 2*besselj(1, pi*u) ./ (pi*u);
end
