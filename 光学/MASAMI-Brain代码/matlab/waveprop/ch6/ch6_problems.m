%% ===== 书 Schmidt PM199 第 6 章（工程 3）=====
%  习题 1-5 复现（原书 6.7 节）
%  与 python/waveprop/ch6/ch6_problems.py 数值行为一致；
%  本机无 MATLAB 时以 Python 实跑结果为准。
%
%  原书习题：
%    1. 高斯激光束的角谱传播：束腰 w0=1mm（R=inf），z2=4m；lambda=1um,
%       N=512, 源面 1cm 网格、观测面 1.5cm 网格。画 y2=0 切片辐照度+相位。
%    2. 圆孔聚焦光束：观测面 = 焦面；lambda=1um, D=1cm, fl=16cm, N=1024,
%       源面 2cm 网格，观测面采样 = 衍射受限光斑直径/100。画焦面辐照度。
%    3. Talbot 成像：tA = (1/2)[1+cos(2*pi*x1/d)]（式 6.93）；观测面取第一
%       Talbot 像面；lambda=1um, d=0.5mm, N=1024, 两面网格 2cm。
%    4. 矩形感兴趣区（宽 Dx, Dy）的模型点源。
%    5. 圆形感兴趣区（直径 D）的模型点源。
%
%  习题 4/5 的推导（式 6.89 + 平移定理，符号按笔记对 (6.92) 的勘误）：
%    矩形窗：U_pt = A (Dx Dy)/(lambda Dz)^2 e^{-i k |r1-rc|^2/(2Dz)}
%            * sinc[Dx(x1-xc)/(lambda Dz)] sinc[Dy(y1-yc)/(lambda Dz)]
%    圆形窗（直径 D）：U_pt = A jinc(|r1-rc|; D/(2 lambda Dz))
%                      e^{-i k |r1-rc|^2/(2Dz)}
%    （两者都是"居中模型平移到 rc"——Fresnel 核平移不变。）
%
%  Python 实跑结果（本脚本按相同参数复现）：
%    习题 1：w 拟合 1.6190 mm（解析 1.6190，2e-12）；R 拟合 6.4674 m（6e-13）；
%            常数相位 -0.9050 = -Gouy；轴上峰值 (w0/w)^2 = 0.381514 精确一致。
%    习题 2：焦面主瓣最小值正落解析第一零点 19.52 um，I(零点) = 3e-5；
%            廓线 vs Airy 最大偏差 1.98e-2（N=1024，陡肩处）-> 5.3e-3（N=2048）。
%    习题 3：zT = 2d^2/lambda = 0.5 m；中央 10 周期最大偏差 7.1e-15。
%    习题 4/5：质心正落 rc；内部 80% 区域 |比值-1| max 7.2e-2 / 6.8e-2。
%
%  MATLAB/Python 差异：
%    - MATLAB 1-based / numpy 0-based（居中样本位置）；
%    - sinc 用本地定义（避免 Signal Processing Toolbox 依赖）；
%    - polyfit/log 拟合与 numpy.polyfit 一致（均按幂次降序返回系数）。
clear; clc; close all;

wvl = 1e-6;
k = 2*pi/wvl;

%% ================= 习题 1：高斯束腰 -> 4 m =================
fprintf('===== 习题 1：高斯束（w0=1 mm 束腰，z=4 m，ang_spec） =====\n');
N1 = 512;  L1 = 1e-2;  L2 = 1.5e-2;
d1 = L1/N1;  d2 = L2/N1;
Dz1 = 4.0;
w0 = 1e-3;
zR = pi*w0^2/wvl;
w_z = w0*sqrt(1 + (Dz1/zR)^2);
R_z = Dz1*(1 + (zR/Dz1)^2);
gouy = atan(Dz1/zR);
n1 = -N1/2 : N1/2-1;
[x1, y1] = meshgrid(n1*d1, n1*d1);
U_in = exp(-(x1.^2 + y1.^2)/w0^2);              % 束腰：场 e^{-r^2/w0^2}
[x2, y2, U_out] = ang_spec_prop(U_in, wvl, d1, d2, Dz1);
I_out = abs(U_out).^2;
I_an = (w0/w_z)^2*exp(-2*(x2.^2 + y2.^2)/w_z^2);
sl = N1/2 + 1;
mfit = (x2(1, :) > 0.05e-3) & (x2(1, :) < 3e-3) & (I_out(sl, :) > 1e-6*max(I_out(:)));
p = polyfit(x2(1, mfit).^2, log(I_out(sl, mfit)), 1);
w_fit = sqrt(-2/p(1));
mph = I_out(sl, :) > 0.05*max(I_out(:));
phi = angle(U_out(sl, mph));
p2 = polyfit(x2(1, mph).^2, phi, 1);
R_fit = k/(2*p2(1));
gouy_fit = mean(phi - p2(1)*x2(1, mph).^2);
fprintf('    zR = %.4f m；解析 w(4m) = %.4f mm，R(4m) = %.4f m，Gouy = %.4f rad\n', ...
    zR, w_z*1e3, R_z, gouy);
fprintf('    m = delta2/delta1 = %.3f\n', d2/d1);
fprintf('    仿真拟合：w = %.4f mm（相对差 %.1e）；R = %.4f m（相对差 %.1e）；常数相位 = %.4f（解析 -Gouy = %.4f）\n', ...
    w_fit*1e3, abs(w_fit - w_z)/w_z, R_fit, abs(R_fit - R_z)/R_z, gouy_fit, -gouy);
fprintf('    轴上峰值辐照度：仿真 %.6f vs 解析 (w0/w)^2 = %.6f\n', ...
    I_out(sl, sl), (w0/w_z)^2);

%% ================= 习题 2：圆孔聚焦 -> 焦面 Airy =================
fprintf('===== 习题 2：圆孔聚焦光束（D=1 cm，fl=16 cm，焦面） =====\n');
N2 = 1024;  L1b = 2e-2;
d1b = L1b/N2;
D_ap = 1e-2;  fl = 0.16;
airy_d = 2.44*wvl*fl/D_ap;
d2b = airy_d/100;
n2 = -N2/2 : N2/2-1;
[x1b, y1b] = meshgrid(n2*d1b, n2*d1b);
ap = circ_(x1b, y1b, D_ap/2);
U_in2 = ap*exp(-1i*k*(x1b.^2 + y1b.^2)/(2*fl));
[x2b, y2b, U_foc] = ang_spec_prop(U_in2, wvl, d1b, d2b, fl);
I_foc = abs(U_foc).^2;
I_foc = I_foc/max(I_foc(:));
r2b = hypot(x2b, y2b);
sl2 = N2/2 + 1;
r_zero = 1.22*wvl*fl/D_ap;
jz = round(r_zero/d2b);
seg = I_foc(sl2, sl2:sl2+59);
jmin = sl2 - 1 + np_argmin(seg);
v = pi*D_ap*abs(x2b(1, :))/(wvl*fl);
airy_sl = airy_fn(v);
mm = v <= 3.8317;
dev = max(abs(I_foc(sl2, mm) - airy_sl(mm)));
fprintf('    delta2 = Airy 直径/100 = %.4f um（m = %.4f）\n', d2b*1e6, d2b/d1b);
fprintf('    Airy 直径 = %.2f um；第一零点 = %.2f um\n', airy_d*1e6, r_zero*1e6);
fprintf('    仿真主瓣最小值位置 = %.2f um；I(解析零点) = %.2e\n', ...
    abs(x2b(1, jmin))*1e6, I_foc(sl2, sl2-1+jz));
fprintf('    廓线 vs Airy（|r|<=零点）最大偏差 = %.2e\n', dev);
% N 翻倍收敛
N2f = 2*N2;  d1f = L1b/N2f;
n2f = -N2f/2 : N2f/2-1;
[x1f, y1f] = meshgrid(n2f*d1f, n2f*d1f);
U_in_f = circ_(x1f, y1f, D_ap/2).*exp(-1i*k*(x1f.^2 + y1f.^2)/(2*fl));
[~, ~, U_foc_f] = ang_spec_prop(U_in_f, wvl, d1f, d2b, fl);
row_f = abs(U_foc_f(N2f/2 + 1, :)).^2;
row_f = row_f/max(row_f);
% 对齐：粗网格 j 对应细网格 j + N2f/2 - N2/2
al = N2f/2 - N2/2;
dev_f = max(abs(row_f(al + (1:N2)) - airy_sl(mm)));
fprintf('    N 翻倍（N=%d）：最大偏差 %.2e（收敛确认）\n', N2f, dev_f);

%% ================= 习题 3：Talbot 成像 =================
fprintf('===== 习题 3：Talbot 成像（d=0.5 mm，第一 Talbot 像） =====\n');
N3 = 1024;  L3 = 2e-2;
d13 = L3/N3;
d_gr = 0.5e-3;
zT = 2*d_gr^2/wvl;
[x3, ~] = meshgrid(n2*d13, n2*d13);
tA = 0.5*(1 + cos(2*pi*x3/d_gr));
[x2t, ~, U_tal] = ang_spec_prop(tA, wvl, d13, d13, zT);
I_tal = abs(U_tal).^2;
I_an3 = abs(tA).^2;      % zT 处各衍射级相位 e^{-i2pi n^2} = 1 -> 精确复制
m10 = abs(x2t(1, :)) <= 5*d_gr;
dev3 = max(abs(I_tal(N3/2 + 1, m10) - I_an3(N3/2 + 1, m10)));
fprintf('    zT = 2 d^2/lambda = %.3f m；m = 1\n', zT);
fprintf('    中央 10 周期：仿真 vs 解析最大偏差 = %.2e\n', dev3);

%% ================= 习题 4/5：矩形/圆形窗模型点源 =================
fprintf('===== 习题 4：矩形窗模型点源（Dx=3 mm, Dy=1.6 mm） =====\n');
N4 = 1024;
Dw4 = 3e-3;
arg4 = Dw4/(wvl*Dz1);
d14 = 1/(10*arg4);
d24 = Dw4/100;
rc4 = [40*d14, -24*d14];
[x14, y14] = meshgrid(n2*d14, n2*d14);
[x24, y24] = meshgrid(n2*d24, n2*d24);
Dx = 3e-3;  Dy = 1.6e-3;
pt4 = pt_source_rect_win(x14, y14, Dx, Dy, wvl, Dz1, rc4);
[~, ~, U4] = ang_spec_prop(pt4, wvl, d14, d24, Dz1);
U_t4 = target_win(x24, y24, 'rect', Dx, Dy, wvl, Dz1, rc4);
[cx, cy, dev_a, dev_p] = compare_win(U4, U_t4, x24, y24, rc4, Dx/2, Dy/2, false);
fprintf('    rc = (%.2f, %.2f) mm；m = %.2f\n', rc4(1)*1e3, rc4(2)*1e3, d24/d14);
fprintf('    传播后质心 = (%.3f, %.3f) mm；内部 80%%：|比值-1| max = %.2e，相位残差 = %.2e rad\n', ...
    cx*1e3, cy*1e3, dev_a, dev_p);
fprintf('===== 习题 5：圆形窗模型点源（D=3 mm） =====\n');
pt5 = pt_source_circ_win(x14, y14, Dw4, wvl, Dz1, rc4);
[~, ~, U5] = ang_spec_prop(pt5, wvl, d14, d24, Dz1);
U_t5 = target_win(x24, y24, 'circ', Dw4, 0, wvl, Dz1, rc4);
[cx, cy, dev_a5, dev_p5] = compare_win(U5, U_t5, x24, y24, rc4, Dw4/2, Dw4/2, true);
fprintf('    传播后质心 = (%.3f, %.3f) mm；内部 80%%：|比值-1| max = %.2e，相位残差 = %.2e rad\n', ...
    cx*1e3, cy*1e3, dev_a5, dev_p5);
fprintf('    jinc 主瓣第一零点 = 0.61 lambda Dz/D = %.4f mm\n', 0.61*wvl*Dz1/Dw4*1e3);

%% 绘图（对应 Python 版图 6-15）
figure('Position', [40 40 1500 820]);
m_ = abs(x2(1, :)) <= 6e-3;
subplot(2, 3, 1);
plot(x2(1, m_)*1e3, I_an(sl, m_), 'k-', x2(1, m_)*1e3, I_out(sl, m_), 'r.');
title('(a) 习题 1：z=4 m 辐照度切片'); grid on;
xlabel('x2 [mm]'); ylabel('辐照度 [a.u.]');
subplot(2, 3, 2);
ph_an = k*x2(1, mph).^2/(2*R_z) - gouy;
ph_an = ph_an - (mean(ph_an - phi));
plot(x2(1, mph)*1e3, phi, 'r.', x2(1, mph)*1e3, ph_an, 'k-');
title('(b) 习题 1：相位切片'); grid on;
xlabel('x2 [mm]'); ylabel('相位 [rad]');
subplot(2, 3, 3);
m_airy = abs(x2b(1, :)) <= 3.5*r_zero;
plot(x2b(1, m_airy)*1e6, airy_sl(m_airy), 'k-', ...
    x2b(1, m_airy)*1e6, I_foc(sl2, m_airy), 'r.');
title('(c) 习题 2：焦面辐照度 vs Airy'); grid on;
xlabel('x2 [um]'); ylabel('归一化辐照度');
subplot(2, 3, 4);
plot(x2t(1, m10)*1e3, I_an3(N3/2 + 1, m10), 'k-', ...
    x2t(1, m10)*1e3, I_tal(N3/2 + 1, m10), 'r.');
title('(d) 习题 3：第一 Talbot 像面（中央 10 周期）'); grid on;
xlabel('x2 [mm]'); ylabel('辐照度 [a.u.]');
subplot(2, 3, 5);
jy4 = N4/2 + round(rc4(2)/d24);
m_s4 = abs(x24(1, :) - rc4(1)) <= 1.6*Dx/2;
plot((x24(1, m_s4) - rc4(1))*1e3, abs(U_t4(jy4, m_s4))/max(abs(U_t4(:))), 'k-', ...
    (x24(1, m_s4) - rc4(1))*1e3, abs(U4(jy4, m_s4))/max(abs(U_t4(:))), 'r.');
title('(e) 习题 4：矩形窗模型点源（y2=yc 切片）'); grid on;
xlabel('(x2-xc) [mm]'); ylabel('归一化 |U|');
subplot(2, 3, 6);
m_s5 = abs(x24(1, :) - rc4(1)) <= 1.6*Dw4/2;
plot((x24(1, m_s5) - rc4(1))*1e3, abs(U_t5(jy4, m_s5))/max(abs(U_t5(:))), 'k-', ...
    (x24(1, m_s5) - rc4(1))*1e3, abs(U5(jy4, m_s5))/max(abs(U_t5(:))), 'r.');
title('(f) 习题 5：圆形窗模型点源（y2=yc 切片）'); grid on;
xlabel('(x2-xc) [mm]'); ylabel('归一化 |U|');

%% ===== 本地函数 =====
function U = pt_source_rect_win(x1, y1, Dx, Dy, wvl, Dz, rc)
% 矩形感兴趣区（宽 Dx, Dy，中心 rc）的模型点源（勘误修正符号）。
A = wvl*Dz;
k = 2*pi/wvl;
argx = Dx/(wvl*Dz);  argy = Dy/(wvl*Dz);
ph = -1i*k*((x1 - rc(1)).^2 + (y1 - rc(2)).^2)/(2*Dz);
U = A*argx*argy*exp(ph) ...
    .* sinc_(argx*(x1 - rc(1))) .* sinc_(argy*(y1 - rc(2)));
end

function U = pt_source_circ_win(x1, y1, D, wvl, Dz, rc)
% 圆形感兴趣区（直径 D，中心 rc）的模型点源：
%   U_pt = A jinc(|r1-rc|; D/(2 lambda Dz)) e^{-i k |r1-rc|^2/(2 Dz)}
A = wvl*Dz;
k = 2*pi/wvl;
rho = hypot(x1 - rc(1), y1 - rc(2));
R = D/(2*wvl*Dz);
% jinc(rho; R) = R*J1(2 pi R rho)/rho（rho=0 处极限 pi R^2）
jn = zeros(size(rho));
nz = rho ~= 0;
jn(nz) = R*besselj(1, 2*pi*R*rho(nz))./rho(nz);
jn(~nz) = pi*R^2;
U = A*jn*exp(-1i*k*rho.^2/(2*Dz));
end

function U = target_win(x2, y2, shape, Dx, Dy, wvl, Dz, rc)
% 式 (6.86) 的目标观测面场（A = lambda*Dz；rect 或 circ 窗）。
A = wvl*Dz;
k = 2*pi/wvl;
if strcmp(shape, 'rect')
    W = rect_(abs(x2 - rc(1)), Dx) .* rect_(abs(y2 - rc(2)), Dy);
else
    W = circ_(x2 - rc(1), y2 - rc(2), Dx/2);
end
U = A*W*exp(1i*k*((x2 - rc(1)).^2 + (y2 - rc(2)).^2)/(2*Dz))/(1i*wvl*Dz);
end

function [cx, cy, dev_a, dev_p] = compare_win(U_sim, U_tgt, x2, y2, rc, hx, hy, circ)
% 窗内质心 + 内部 80% 区域的复场偏差统计（rect: max-norm；circ: 2-norm）。
I = abs(U_sim).^2;
w = I > 0.05*max(I(:));
cx = sum(x2(w).*I(w))/sum(I(w));
cy = sum(y2(w).*I(w))/sum(I(w));
if circ
    m_in = ((x2 - rc(1))/hx).^2 + ((y2 - rc(1))/hy).^2 <= 0.64;
else
    m_in = (abs(x2 - rc(1)) <= 0.8*hx) & (abs(y2 - rc(2)) <= 0.8*hy);
end
ratio = U_sim(m_in)./U_tgt(m_in);
dev_a = max(abs(abs(ratio) - 1));
dev_p = max(abs(angle(ratio) - mean(angle(ratio))));
end

function i = np_argmin(v)
% numpy.argmin 语义：第一个最小值索引（MATLAB min 的第一个输出同义）
[~, i] = min(v);
end

function a = airy_fn(v)
% Airy 函数 [2 J1(v)/v]^2（v=0 处极限 1）
a = zeros(size(v));
nz = v > 1e-12;
a(nz) = (2*besselj(1, v(nz))./v(nz)).^2;
a(~nz) = 1;
end

function out = rect_(x, width)
% 宽度 width 的矩形窗
out = double(abs(x) < width/2) + 0.5*double(abs(x) == width/2);
end

function out = circ_(x, y, R)
% 半径 R 圆盘（半径参数）
out = double(hypot(x, y) <= R);
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
