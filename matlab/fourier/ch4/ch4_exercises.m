% ===== 书 4.6 节习题 4.2-4.5 参考解答（数值部分，与 ch4_exercises.py 同构） =====
% 4.2 三种衍射解的适用距离（解析估算 + 打印）；
% 4.3 验证菲涅耳传递函数 (4.24) = 菲涅耳脉冲响应 (4.23) 的傅里叶变换；
% 4.4 三个孔径的解析夫琅禾费辐照度：(a) 矩形孔 (b) 环形孔 (c) 双圆孔；
% 4.5 双圆孔两路振幅不同（A1=1, A2=0.4）。
% 解析式统一用书 4.4.3 节写法：场 = e^{jkz}/(j lam z) * e^{jk(x^2+y^2)/(2z)}
%   * FT{U1}(fX = x/(lam z), fY = y/(lam z))；辐照度按式 (4.8) 取模方。
% 依赖：同文件夹 jinc.m（书 4.5 节建）；sinc 用文件末尾的局部函数 sincn
% （基础 MATLAB 没有 sinc——它是信号处理工具箱函数，定义与 np.sinc 相同）。
% MATLAB/Python 差异：MATLAB 索引 1-based，y=0 的剖面在 I2(M/2+1,:)；
%   Python 0-based 在 I2[M//2]。
clear; clc; close all;

lam = 0.633e-6;    % 习题 4.4 起统一用 He-Ne 波长（习题 4.2 另用 0.5 um）

% ================= 习题 4.2：三种衍射解的适用距离 =================
w42   = 0.5e-3;      % 孔径半径（直径 1 mm）
lam42 = 0.5e-6;      % 习题指定波长
% 式 (4.27)：z^3 >> (pi/(4 lam)) * [(x-xi)^2+(y-eta)^2]^2_max
% 观察区取孔径的几何影区（直径 1 mm）-> x-xi、y-eta 各可达 1 mm
Dmax2  = (2*w42)^2 + (2*w42)^2;
z_eq27 = (pi/(4*lam42) * Dmax2^2)^(1/3);
z_nf1  = w42^2/lam42;                       % NF = 1
z_nf01 = 10*w42^2/lam42;                    % NF = 0.1（书 4.5 节工程判据）
z_fh   = pi*(w42^2 + w42^2)/lam42;          % 式 (4.29)：z >> k(xi^2+eta^2)_max/2
fprintf('========== 习题 4.2（直径 1 mm 孔径，lam = 0.5 um） ==========\n');
fprintf('RS：只需 r >> lam -> 任何宏观距离都满足（z > ~1 mm 即可）\n');
fprintf('Fresnel：式(4.27)（观察区取几何影区）z >> %.1f cm；NF 视角：NF<1 对应 z > %.2f m，平滑孔径可到 NF~20-30（z ~ %.3f m）\n', ...
        z_eq27*1e2, z_nf1, z_nf1/25);
fprintf('Fraunhofer：式(4.29) z >> %.2f m；工程口径 NF<0.1 -> z > %.1f m\n', z_fh, z_nf01);

% ================= 习题 4.3：H = FT{h} 的数值验证 =================
% h(x,y) = e^{jkz}/(j lam z) * exp[j k (x^2+y^2)/(2z)]   （式 4.23）
% H(fX,fY) = e^{jkz} * exp[-j pi lam z (fX^2+fY^2)]      （式 4.24）
% 比较频带 |f| <= L/(2 lam z)（驻点带），且要低于 Nyquist（第 5 章判据预览）。
L3 = 0.02; M3 = 512; z3 = 2;
dx3 = L3/M3;
x3 = -L3/2 : dx3 : L3/2-dx3;
f3 = -1/(2*dx3) : 1/L3 : 1/(2*dx3)-1/L3;    % 式 (2.22)
[X3, Y3] = meshgrid(x3, x3);
[FX3, FY3] = meshgrid(f3, f3);
k3 = 2*pi/lam;
h = exp(1j*k3*z3)/(1j*lam*z3) * exp(1j*k3*(X3.^2+Y3.^2)/(2*z3));
H_num = fftshift(fft2(ifftshift(h)))*dx3^2;   % FT ≈ dx^2 · FFT
H_an  = exp(1j*k3*z3) * exp(-1j*pi*lam*z3*(FX3.^2+FY3.^2));

band = L3/(2*lam*z3);
inb8 = (abs(FX3) <= 0.8*band) & (abs(FY3) <= 0.8*band);
inb5 = (abs(FX3) <= 0.5*band) & (abs(FY3) <= 0.5*band);
dphi8 = angle(H_num(inb8).*conj(H_an(inb8)));
dphi5 = angle(H_num(inb5).*conj(H_an(inb5)));
fprintf('\n========== 习题 4.3（L=0.02 m, M=512, z=2 m） ==========\n');
fprintf('驻点带 |f| <= L/(2 lam z) = %.3e cyc/m（Nyquist = %.3e）\n', band, 1/(2*dx3));
fprintf('0.8 频带内最大缠绕相位差 = %.3f rad\n', max(abs(dphi8)));
fprintf('0.5 频带内最大缠绕相位差 = %.3f rad（残差来自常数幅度 chirp 的窗口截断）\n', ...
        max(abs(dphi5)));

figure('Position', [60 60 1050 420]);
subplot(1,2,1);
sel = abs(f3) <= 0.8*band;
dphi_line = angle(H_num(M3/2+1,:).*conj(H_an(M3/2+1,:)));
plot(f3(sel), dphi_line(sel));
xlabel('f_X (cyc/m)'); ylabel('wrapped phase diff (rad)');
title('Ex 4.3: \angle(H_{num}/H_{an}) along f_X (f_Y=0)');
grid on;
subplot(1,2,2);
sel2 = abs(f3) <= 400;
plot(f3(sel2), unwrap(angle(H_num(M3/2+1,sel2)))); hold on;
plot(f3(sel2), unwrap(angle(H_an(M3/2+1,sel2))), ':');
xlabel('f_X (cyc/m)'); ylabel('unwrapped phase (rad)');
title('low-freq window: the two chirps coincide');
legend('FFT\{h\}', 'H Eq.(4.24)', 'FontSize', 8); grid on;
saveas(gcf, 'out_ex43.png');

% ================= 习题 4.4：三个孔径的解析夫琅禾费图样 =================
% (a) 矩形孔 wX=0.1 mm, wY=0.05 mm, z=5 m,  L=0.16 m
% (b) 环形孔 w1=1 mm, w2=0.2 mm,             z=50 m,  L=0.2 m
% (c) 双圆孔 w=1 mm, 中心距 s=4 mm,          z=500 m, L=0.6 m
% z 都按式 (4.29) 取 ~10 倍以上裕量；L ≈ 5 倍中心特征宽度。
fprintf('\n========== 习题 4.4 ==========\n');

% (a) 矩形孔：I2 = (4wXwY/(lam z))^2 sinc^2(2wX x/(lam z)) sinc^2(2wY y/(lam z))
wX = 0.1e-3; wY = 0.05e-3; za = 5; La = 0.16; Ma = 512;
dxa = La/Ma; xa = -La/2 : dxa : La/2-dxa;
[Xa, Ya] = meshgrid(xa, xa);
lza = lam*za;
Ia = (4*wX*wY/lza)^2 .* sincn(2*wX*Xa/lza).^2 .* sincn(2*wY*Ya/lza).^2;
fprintf('(a) z=5 m：式(4.29)要求 z >> %.3f m（裕量 %.0f 倍）；峰 = %.3e，实测 %.3e\n', ...
        pi*(wX^2+wY^2)/lam, 5/(pi*(wX^2+wY^2)/lam), (4*wX*wY/lza)^2, ...
        Ia(Ma/2+1, Ma/2+1));
fprintf('    sinc 第一零点 x = lam z/(2wX) = %.4f m（y 向 %.4f m）\n', ...
        lza/(2*wX), lza/(2*wY));

% (b) 环形孔：amp = (w1^2/(lam z)) [jinc(w1 rho/(lam z)) - (w2/w1)^2 jinc(w2 rho/(lam z))]
w1 = 1e-3; w2 = 0.2e-3; zb = 50; Lb = 0.2; Mb = 512;
dxb = Lb/Mb; xb = -Lb/2 : dxb : Lb/2-dxb;
[Xb, Yb] = meshgrid(xb, xb);
lzb = lam*zb;
rho_b = sqrt(Xb.^2+Yb.^2);
amp_b = w1^2/lzb * (jinc(w1*rho_b/lzb) - (w2/w1)^2 * jinc(w2*rho_b/lzb));
Ib = amp_b.^2;
fprintf('(b) z=50 m：式(4.29)要求 z >> %.2f m（裕量 %.0f 倍）；NF = %.4f\n', ...
        pi*(w1^2+w2^2)/lam, 50/(pi*(w1^2+w2^2)/lam), w1^2/lzb);
fprintf('    峰 = (pi(w1^2-w2^2)/(lam z))^2 = %.3e，实测 %.3e（无遮拦 %.3e）\n', ...
        (pi*(w1^2-w2^2)/lzb)^2, Ib(Mb/2+1,Mb/2+1), (pi*w1^2/lzb)^2);
fprintf('    外环第一零点 ≈ 0.61*lam*z/w1 = %.5f m\n', 0.61*lzb/w1);

% (c) 双圆孔：I2 = (2 w^2/(lam z))^2 cos^2(pi s x/(lam z)) jinc^2(w rho/(lam z))
wc = 1e-3; sc = 4e-3; zc = 500; Lc = 0.6; Mc = 512;
dxc = Lc/Mc; xc = -Lc/2 : dxc : Lc/2-dxc;
[Xc, Yc] = meshgrid(xc, xc);
lzc = lam*zc;
rho_c = sqrt(Xc.^2+Yc.^2);
env_c = wc^2/lzc * jinc(wc*rho_c/lzc);        % 单孔包络（含 1/(lam z) 前置）
Ic = (2*env_c).^2 .* cos(pi*sc*Xc/lzc).^2;
fprintf('(c) z=500 m：式(4.29)要求 z >> %.1f m（裕量 %.0f 倍）；NF = %.4f\n', ...
        pi*((sc/2+wc)^2+wc^2)/lam, 500/(pi*((sc/2+wc)^2+wc^2)/lam), wc^2/lzc);
fprintf('    峰 = (2 pi w^2/(lam z))^2 = %.3e，实测 %.3e\n', ...
        (2*pi*wc^2/lzc)^2, Ic(Mc/2+1,Mc/2+1));
fprintf('    包络第一零点 0.61*lam*z/w = %.4f m；条纹周期 lam*z/s = %.5f m\n', ...
        0.61*lzc/wc, lzc/sc);

figure('Position', [60 60 1250 720]);
subplot(2,3,1);
imagesc(xa, xa, nthroot(Ia,3)); colormap('gray'); axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)'); title('(a) rect: I_2^{1/3} image');
subplot(2,3,2);
plot(xa, Ia(Ma/2+1,:)); grid on;
xlabel('x (m)'); ylabel('Irradiance');
title('(a) profile: sinc^2 \times sinc^2');
subplot(2,3,3);
plot(xb, Ib(Mb/2+1,:)); hold on;
plot(xb, (w1^2/lzb)^2 * jinc(w1*xb/lzb).^2, ':'); grid on;
xlabel('x (m)'); ylabel('Irradiance'); title('(b) profile: annulus vs unobscured');
legend('annulus', 'unobscured disk', 'FontSize', 8);
subplot(2,3,4);
imagesc(xb, xb, nthroot(Ib,3)); colormap('gray'); axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)'); title('(b) annulus: I_2^{1/3} image');
subplot(2,3,5);
imagesc(xc, xc, nthroot(Ic,3)); colormap('gray'); axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)'); title('(c) pair: I_2^{1/3} image');
subplot(2,3,6);
plot(xc, Ic(Mc/2+1,:)); grid on;
xlabel('x (m)'); ylabel('Irradiance');
title('(c) profile: fringes \times Airy envelope');
saveas(gcf, 'out_ex44.png');

% ================= 习题 4.5：两路振幅不同（A1=1, A2=0.4） =================
% I2 = (w^2/(lam z))^2 jinc^2(...) * (A1^2 + A2^2 + 2 A1 A2 cos(2 pi s x/(lam z)))
% （书提示：先取模方，再把两个复指数合并成余弦）
A1 = 1; A2 = 0.4;
Ie = (wc^2/lzc)^2 .* jinc(wc*rho_c/lzc).^2 .* ...
     (A1^2 + A2^2 + 2*A1*A2*cos(2*pi*sc*Xc/lzc));
fprintf('\n========== 习题 4.5（A1=1, A2=0.4, z=500 m） ==========\n');
fprintf('峰值 = (w^2/(lam z))^2 pi^2 (A1+A2)^2 = %.3e，实测 %.3e\n', ...
        (wc^2/lzc)^2*pi^2*(A1+A2)^2, max(Ie(Mc/2+1,:)));
fprintf('条纹可见度 V = 2A1A2/(A1^2+A2^2) = %.4f；Imax/Imin = ((A1+A2)/(A1-A2))^2 = %.2f\n', ...
        2*A1*A2/(A1^2+A2^2), ((A1+A2)/(A1-A2))^2);

figure('Position', [80 80 1050 420]);
subplot(1,2,1);
imagesc(xc, xc, nthroot(Ie,3)); colormap('gray'); axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title('Ex 4.5: A_1=1, A_2=0.4 — I_2^{1/3} image');
subplot(1,2,2);
plot(xc, Ie(Mc/2+1,:)); grid on;
xlabel('x (m)'); ylabel('Irradiance');
title('profile: visibility 2A_1A_2/(A_1^2+A_2^2) = 0.69');
saveas(gcf, 'out_ex45.png');
fprintf('\n已保存 out_ex43.png, out_ex44.png, out_ex45.png\n');

% ---------- 局部函数 ----------
function s = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)，u=0 处极限 1（表 1.2；np.sinc 同定义）
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end
