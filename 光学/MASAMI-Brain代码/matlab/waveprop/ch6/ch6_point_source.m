%% ===== 书 Schmidt PM199 第 6 章（工程 2）=====
%  点源模型：sinc 模型与 sinc-Gaussian 模型（原书 6.6 节，Listings 6.8-6.9）
%  与 python/waveprop/ch6/ch6_point_source.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准。
%
%  原书要点：
%    真点源 = Dirac delta（式 6.82），谱无限宽 —— 离散网格必须用带限模型：
%      Martin-Flatté（式 6.83）：exp(-r^2/2sigma^2) exp(-i r^2/2x0^2)；
%      Flatté et al.（式 6.84）：exp(-r^2/2sigma^2) cos^2(r^2/2rho^2)；
%      本书做法（逆向设计）：要求观测面给出"窗口加权的抛物面相位球面波"
%      （式 6.85/6.86），反解源面模型（式 6.87-6.92）：
%        方窗（宽 D，中心 r_c）：sinc 模型（式 6.92，Listing 6.8）；
%        乘 Gaussian 抑制旁瓣 -> sinc-Gaussian 模型（Listing 6.9）。
%
%  原书勘误（详见笔记）：
%    式 (6.92) 的两个指数因子符号印反：e^{+i k rc^2/(2Dz)} e^{-i(k/Dz) rc.r1}
%    应为 e^{-i k rc^2/(2Dz)} e^{+i(k/Dz) rc.r1}。修正后模型恰为居中模型的
%    平移（Fresnel 核平移不变）；书印版本的额外相位斜坡等价于倾斜入射，
%    把窗口推到 -rc —— Python 实测质心 (1.20,0.80) vs (-1.20,-0.80) mm。
%
%  验证（编号与 Python 版一致）：
%   (0) 源面模型峰值辐照度 = A^2 arg^4 = 4.096 kW/m^2（原书图 6.7 ~4.1）；
%   (1) Listing 6.8 复现：观测面窗内平顶 ~1 W/m^2（图 6.8）、相位混叠（图 6.9）；
%   (2) Listing 6.9 复现：sinc-Gaussian 的相位混叠减轻（图 6.10-6.12）；
%   (3) Parseval 能量守恒；
%   (4) 勘误 (6.92) 对照：离轴窗口质心 +rc（修正）vs -rc（书印）。
%
%  MATLAB/Python 差异：
%    - MATLAB 1-based / numpy 0-based；sinc 在 Signal Processing Toolbox
%      （sin(pi x)/(pi x)，与 numpy 一致）；本脚本自带定义避免工具箱依赖。
clear; clc; close all;

%% 原书 Listing 6.8/6.9 的参数
D   = 8e-3;          % diameter of the observation aperture [m]
wvl = 1e-6;          % optical wavelength [m]
Dz  = 1;             % propagation distance [m]
arg = D/(wvl*Dz);
d1  = 1/(10*arg);    % source-plane grid spacing（中心瓣 10 采样/半瓣）
d2  = D/100;         % observation-plane grid spacing
N   = 1024;
n  = -N/2 : N/2-1;
[x1, y1] = meshgrid(n*d1, n*d1);
[x2, y2] = meshgrid(n*d2, n*d2);
r2sq = x2.^2 + y2.^2;

%% 验证 0：源面模型峰值辐照度
pt  = pt_source_sinc(x1, y1, D, wvl, Dz);
ptg = pt_source_sinc_gauss(x1, y1, D, wvl, Dz);
I_pt = abs(pt).^2;  I_ptg = abs(ptg).^2;
fprintf('===== 验证 0：源面模型峰值辐照度 =====\n');
fprintf('    A = lambda*Dz = %.2e m；arg = %.0f /m；delta1 = %.2f um，delta2 = %.0f um\n', ...
    wvl*Dz, arg, d1*1e6, d2*1e6);
fprintf('    sinc 模型：源面峰值辐照度 = %.3f kW/m^2（解析 A^2 arg^4 = %.3f）\n', ...
    max(I_pt(:))/1e3, (wvl*Dz)^2*arg^4/1e3);
fprintf('    sinc-Gauss：源面峰值辐照度 = %.3f kW/m^2（高斯因子在 r=0 处为 1，压旁瓣不压峰值）\n', ...
    max(I_ptg(:))/1e3);

%% 验证 1/2：传播与观测面统计（原书图 6.8/6.9 与 6.11/6.12）
[~, ~, U_ob]  = ang_spec_prop(pt,  wvl, d1, d2, Dz);
[~, ~, U_obg] = ang_spec_prop(ptg, wvl, d1, d2, Dz);
I_ob = abs(U_ob).^2;  I_obg = abs(U_obg).^2;
E_in  = sum(abs(pt(:)).^2)*d1^2;
E_out = sum(I_ob(:))*d2^2;
fprintf('===== 验证 1/2：传播与观测面统计 =====\n');
fprintf('    Parseval：sum|U|^2 d2^2 / sum|pt|^2 d1^2 = %.6f\n', E_out/E_in);
m_win = r2sq <= (0.95*D/2)^2;
fprintf('    sinc：窗内 mean = %.4f W/m^2，std/mean = %.2e\n', ...
    mean(I_ob(m_win)), std(I_ob(m_win))/mean(I_ob(m_win)));
fprintf('    sinc-Gauss：窗内 mean = %.4f W/m^2，std/mean = %.2e\n', ...
    mean(I_obg(m_win)), std(I_obg(m_win))/mean(I_obg(m_win)));
Uan = exp(1i*2*pi/wvl/(2*Dz)*r2sq)/(1i*wvl*Dz);   % 抛物线相位目标（无窗）
[r_in_i, r_out_i] = deal(0.95*D/2, 1.10*D/2);
dph  = angle(U_ob .* conj(Uan));
dphg = angle(U_obg .* conj(Uan));
dph  = dph  - mean(dph(r2sq <= r_in_i^2));
dphg = dphg - mean(dphg(r2sq <= r_in_i^2));
m_out = (r2sq >= r_out_i^2) & (r2sq <= 6e-3^2);
fprintf('    sinc：相位误差 RMS 窗内 = %.4f rad；窗外环带 = %.4f rad（混叠）；环带峰值 |U| = %.3f\n', ...
    sqrt(mean(dph(r2sq <= r_in_i^2).^2)), sqrt(mean(dph(m_out).^2)), ...
    max(abs(U_ob(m_out))));
fprintf('    sinc-Gauss：相位误差 RMS 窗内 = %.4f rad；窗外环带 = %.4f rad；环带峰值 |U| = %.3f\n', ...
    sqrt(mean(dphg(r2sq <= r_in_i^2).^2)), sqrt(mean(dphg(m_out).^2)), ...
    max(abs(U_obg(m_out))));

%% 验证 3：勘误 (6.92)：离轴窗口
Dw   = 4e-3;
argw = Dw/(wvl*Dz);
d1w  = 1/(10*argw);
rc   = [48*d1w, 32*d1w];          % = (1.2 mm, 0.8 mm)，delta1 整数倍
[x1w, y1w] = meshgrid(n*d1w, n*d1w);
[x2w, y2w] = meshgrid(n*d2, n*d2);
pt_ok = pt_source_offaxis(x1w, y1w, Dw, wvl, Dz, rc, true);
pt_bad = pt_source_offaxis(x1w, y1w, Dw, wvl, Dz, rc, false);
[~, ~, U_ok]  = ang_spec_prop(pt_ok,  wvl, d1w, d2, Dz);
[~, ~, U_bad] = ang_spec_prop(pt_bad, wvl, d1w, d2, Dz);
I_ok = abs(U_ok).^2;  I_bad = abs(U_bad).^2;
[c_ok, c_bad] = deal(centroid_of(I_ok, x2w, y2w), centroid_of(I_bad, x2w, y2w));
fprintf('===== 验证 3：勘误 (6.92)：离轴窗口 rc = (%.2f, %.2f) mm =====\n', rc(1)*1e3, rc(2)*1e3);
fprintf('    修正模型：传播后质心 = (%.2f, %.2f) mm —— 落在 +rc\n', c_ok(1)*1e3, c_ok(2)*1e3);
fprintf('    书印模型：传播后质心 = (%.2f, %.2f) mm —— 落在 -rc（相位斜坡 = 倾斜入射）\n', ...
    c_bad(1)*1e3, c_bad(2)*1e3);
U_tgt = target_windowed_wave(x2w, y2w, Dw, wvl, Dz, rc);
m_t = abs(U_tgt) > 0.5*max(abs(U_tgt(:)));
m_in = ((x2w - rc(1))/(Dw/2)).^2 + ((y2w - rc(2))/(Dw/2)).^2 <= 0.64;
ratio = U_ok(m_t & m_in)./U_tgt(m_t & m_in);
fprintf('    修正模型 vs 目标场（内部 80%%）：|比值-1| max = %.2e；相位残差 max = %.2e rad\n', ...
    max(abs(abs(ratio) - 1)), max(abs(angle(ratio) - mean(angle(ratio)))));

%% 绘图（对应 Python 版图 6-14）
sl = N/2 + 1;
figure('Position', [40 40 1500 820]);
subplot(2, 3, 1);
m1 = abs(x1(1, :)) <= 0.75e-3;
plot(x1(1, m1)*1e3, I_pt(sl, m1)/1e3, 'k-');
title('(a) sinc 模型点源（Listing 6.8）'); grid on;
xlabel('x1 [mm]'); ylabel('辐照度 [kW/m^2]');
subplot(2, 3, 2);
plot(x1(1, m1)*1e3, I_ptg(sl, m1)/1e3, 'k-');
title('(b) sinc-Gaussian 模型（Listing 6.9）'); grid on;
xlabel('x1 [mm]'); ylabel('辐照度 [kW/m^2]');
subplot(2, 3, 3);
m2 = abs(x2(1, :)) <= 8e-3;
plot(x2(1, m2)*1e3, I_ob(sl, m2), 'b-', x2(1, m2)*1e3, I_obg(sl, m2), 'r-.');
xline(D/2*1e3, ':'); xline(-D/2*1e3, ':');
title('(c) 观测面辐照度切片：窗内平顶 ~1 W/m^2'); grid on;
xlabel('x2 [mm]'); ylabel('辐照度 [W/m^2]');
subplot(2, 3, 4);
plot(x2(1, m2)*1e3, angle(U_ob(sl, m2)), 'b.', x2(1, m2)*1e3, ...
    angle(exp(1i*(2*pi/wvl/(2*Dz)*x2(1, m2).^2 - pi/2))), 'k-');
title('(d) sinc 模型观测面相位（窗外混叠）'); grid on;
xlabel('x2 [mm]'); ylabel('相位 [rad]');
subplot(2, 3, 5);
plot(x2(1, m2)*1e3, angle(U_obg(sl, m2)), 'r.', x2(1, m2)*1e3, ...
    angle(exp(1i*(2*pi/wvl/(2*Dz)*x2(1, m2).^2 - pi/2))), 'k-');
title('(e) sinc-Gaussian 模型观测面相位（混叠减轻）'); grid on;
xlabel('x2 [mm]'); ylabel('相位 [rad]');
subplot(2, 3, 6);
m4 = abs(x2(1, :)) <= 6e-3;
jy = N/2 + round(rc(2)/d2);
plot((x2w(1, m4) - rc(1))*1e3, I_ok(jy, m4)/max(I_ok(:)), 'b-', ...
    (x2w(1, m4) - rc(1))*1e3, I_bad(jy, m4)/max(I_ok(:)), 'r--');
xline(rc(1)*1e3, 'b:'); xline(-rc(1)*1e3, 'r:');
title('(f) 勘误 (6.92)：书印模型窗口推到 -rc'); grid on;
xlabel('(x2-xc) [mm]'); ylabel('归一化辐照度');

%% ===== 本地函数 =====
function U = pt_source_sinc(x1, y1, D, wvl, Dz)
% sinc 模型点源（式 6.92，rc = 0；Listing 6.8）。A = lambda*Dz 定标窗内振幅 1。
A = wvl*Dz;
arg = D/(wvl*Dz);
r1sq = x1.^2 + y1.^2;
U = A*arg^2*exp(-1i*2*pi/wvl/(2*Dz)*r1sq) ...
    .* sinc_(arg*x1) .* sinc_(arg*y1);
end

function U = pt_source_sinc_gauss(x1, y1, D, wvl, Dz)
% sinc-Gaussian 模型点源（Listing 6.9：sinc 模型乘 exp(-(arg/4 r1)^2)）。
A = wvl*Dz;
arg = D/(wvl*Dz);
r1sq = x1.^2 + y1.^2;
U = A*arg^2*exp(-1i*2*pi/wvl/(2*Dz)*r1sq) ...
    .* sinc_(arg*x1) .* sinc_(arg*y1) ...
    .* exp(-(arg/4*sqrt(r1sq)).^2);
end

function U = pt_source_offaxis(x1, y1, D, wvl, Dz, rc, corrected)
% 式 (6.92) 的离轴推广（窗口中心 rc，方窗宽 D）。
% corrected=true（勘误修正版）：恰为居中模型的平移：
%   U_pt = A arg^2 e^{-i k |r1-rc|^2/(2Dz)} sinc[arg(x1-xc)] sinc[arg(y1-yc)]
% corrected=false（书印 6.92）：指数符号印反。
A = wvl*Dz;
arg = D/(wvl*Dz);
k = 2*pi/wvl;
xc = rc(1);  yc = rc(2);
sinc_env = sinc_(arg*(x1 - xc)).*sinc_(arg*(y1 - yc));
if corrected
    ph = -1i*k*((x1 - xc).^2 + (y1 - yc).^2)/(2*Dz);
else
    ph = -1i*k*(x1.^2 + y1.^2)/(2*Dz) ...
        + 1i*k*(xc^2 + yc^2)/(2*Dz) ...
        - 1i*k/Dz*(xc*x1 + yc*y1);
end
U = A*arg^2*exp(ph).*sinc_env;
end

function U = target_windowed_wave(x2, y2, D, wvl, Dz, rc)
% 目标观测面场（式 6.86/6.90，省略 e^{ikDz}）：A = lambda*Dz，方窗。
A = wvl*Dz;
k = 2*pi/wvl;
W = rect_(abs(x2 - rc(1)), D) .* rect_(abs(y2 - rc(2)), D);
U = A*W*exp(1i*k*((x2 - rc(1)).^2 + (y2 - rc(2)).^2)/(2*Dz))/(1i*wvl*Dz);
end

function c = centroid_of(I, x2, y2)
% 显著区域（|U|^2 > 5% 峰值）的辐照度质心
w = I > 0.05*max(I(:));
c = [sum(x2(w).*I(w))/sum(I(w)), sum(y2(w).*I(w))/sum(I(w))];
end

function out = rect_(x, width)
% 宽度 width 的矩形窗
out = double(abs(x) < width/2) + 0.5*double(abs(x) == width/2);
end

function out = sinc_(x)
% 归一化 sinc（避免依赖 Signal Processing Toolbox）
out = ones(size(x));
nz = x ~= 0;
out(nz) = sin(pi*x(nz))./(pi*x(nz));
end

function G = ft2(g, delta)
% 原书 Listing 2.5：定标正确的二维 DFT（乘 delta^2）
G = fftshift(fft2(fftshift(g)))*delta^2;
end

function g = ift2(G, delta_f)
% 原书 Listing 2.6：定标正确的二维逆 DFT（乘 (N*delta_f)^2）
N = size(G, 1);
g = fftshift(ifft2(fftshift(G)))*(N*delta_f)^2;
end

function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
% 书 Listing 6.5：角谱传播（式 6.50 的 m 分支）；与 ch6_fresnel_prop.m 相同
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
