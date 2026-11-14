% ===== 书 5.7 节习题 5.1 / 5.3 / 5.5 / 5.6 的数值解答(与 ch5_exercises.py 同构) =====
% 5.1 圆孔 + 临界采样参数设计(L1 = sqrt(M*lambda*z));
% 5.3 功率守恒诊断:P = sum(sum(I))*dx^2,源面 vs 观察面;
% 5.5 高斯光束:z_R = pi*w0^2/lambda,临界距离 z = L^2/(M*lambda),
%     数值传播 vs 式 (5.26) 解析辐照度;
% 5.6 split-step 等价性(式 (5.15)):10 x TF(2000 m) vs 单次 TF(20000 m)。
% 依赖:同文件夹 propTF.m / propIR.m(书 5.1/5.2 节)。
% MATLAB/Python 差异:MATLAB 索引 1-based(剖面在 (M/2+1,:));
%   MATLAB R2020b 起有 (2)^(1/2) 写法;幂用 ^,逐点用 .^。
clear; clc;

% ================= 习题 5.1:圆孔 + 临界采样 =================
% w = 0.05 m, lambda = 0.5 um, z = 1000 m, M = 500,临界采样。
w1 = 0.05; lam1 = 0.5e-6; z1c = 1000; M1 = 500;
L1c = sqrt(M1*lam1*z1c);      % 临界:dx = L1/M = lambda z/L1 -> L1^2 = M*lambda*z
dx1c = L1c/M1;
FN = 1/(2*dx1c);              % Nyquist 频率
B1 = 5.98/w1;                 % 有效带宽估计(式 5.10)
NF1 = w1^2/(lam1*z1c);        % 菲涅耳数
fprintf('========== 习题 5.1(圆孔 w=0.05 m, lambda=0.5 um, z=1000 m, M=500) ==========\n');
fprintf('(a) L1 = sqrt(M*lambda*z) = %.3f m,dx = %.2f mm,FN = %.0f cyc/m\n', ...
        L1c, dx1c*1e3, FN);
fprintf('(b) B1 = 5.98/w = %.1f < FN = %.0f ✓;跨直径样本数 = %.0f\n', ...
        B1, FN, 2*w1/dx1c);
fprintf('(c) NF = %.1f(z_NF=1 为 %.0f m)-> 处于菲涅耳区\n', NF1, w1^2/lam1);

L1 = L1c;
x1 = -L1/2:dx1c:L1/2-dx1c;
[X1, Y1] = meshgrid(x1, x1);
u1c = double(hypot(X1, Y1) < w1);   % 圆孔(半径 = 分母的约定,见 ch4 勘误 1)
I1c = abs(u1c).^2;
fprintf('(d) 圆孔功率 = pi*w^2 = %.6f,离散 %.6f\n', pi*w1^2, sum(sum(I1c))*dx1c^2);

figure('Position', [60 60 1250 380]);
for i = 1:3
    zt = [500 1000 2000];
    I2 = abs(propTF(u1c, L1, lam1, zt(i))).^2;
    subplot(1,3,i);
    imagesc(x1, x1, nthroot(I2,3)); colormap('gray'); axis square; axis xy;
    xlabel('x (m)'); ylabel('y (m)');
    ratio = lam1*zt(i)/(L1*dx1c);   % lambda z/(L dx):=1 临界
    if ratio < 0.99
        reg = 'TF over';
    elseif ratio > 1.01
        reg = 'TF under';
    else
        reg = 'critical';
    end
    title(sprintf('z = %d m, %s', zt(i), reg));
end
saveas(gcf, 'out_ex5_1.png');

% ================= 习题 5.3:功率守恒诊断 =================
% P = ∫∫ I dxdy ≈ sum(sum(I))*dx^2(无吸收/散射时应守恒;书 5.3 题)。
Lq = 0.5; Mq = 250; lamq = 0.5e-6; wq = 0.051;
dxq = Lq/Mq;
xq = -Lq/2:dxq:Lq/2-dxq;
[Xq, Yq] = meshgrid(xq, xq);
u1q = double(abs(Xq) < wq & abs(Yq) < wq);
P1 = sum(sum(abs(u1q).^2))*dxq^2;
fprintf('\n========== 习题 5.3(功率守恒,源面功率 = %.6f) ==========\n', P1);
fprintf('   z (m) |     P_TF |   P_TF/P1 |     P_IR |   P_IR/P1\n');
for z = [1000 2000 4000 20000]
    u2tf = propTF(u1q, Lq, lamq, z);
    u2ir = propIR(u1q, Lq, lamq, z);
    Ptf = sum(sum(abs(u2tf).^2))*dxq^2;
    Pir = sum(sum(abs(u2ir).^2))*dxq^2;
    fprintf('%8d | %8.6f | %9.5f | %8.6f | %9.5f\n', z, Ptf, Ptf/P1, Pir, Pir/P1);
end
fprintf('(TF 全程守恒;IR@1000 的周期副本把能量复制 4 份 -> P/P1 = 4;\n');
fprintf(' IR@20000 压掉 L/(2 lambda z) 外的源频谱 -> 功率略降)\n');

% ================= 习题 5.5:高斯光束 =================
w0 = 1e-3; lam5 = 0.633e-6; A0 = 1.0;   % 式 (5.29) 参数
L5 = 0.015; M5 = 250;                   % 边长 15 mm,250x250
dx5 = L5/M5;
zR = pi*w0^2/lam5;                      % 式 (5.28):Rayleigh 范围
zcrit = L5^2/(M5*lam5);                 % 临界采样距离(式 5.14)
x5 = -L5/2:dx5:L5/2-dx5;
[X5, Y5] = meshgrid(x5, x5);
U0 = A0*exp(-(X5.^2+Y5.^2)/w0^2);       % 式 (5.29)
fprintf('\n========== 习题 5.5(高斯光束 w0=1 mm, lambda=0.633 um, L=15 mm, M=250) ==========\n');
fprintf('z_R = pi*w0^2/lambda = %.3f m;临界采样距离 z = L^2/(M lambda) = %.3f m\n', ...
        zR, zcrit);
avail = L5/(2*lam5*10);
fprintf('(b) z=10 m 可用带宽 L/(2 lambda z) = %.0f cyc/m(高斯谱带宽远小于此,判据满足)\n', ...
        avail);
figure('Position', [60 60 1250 380]);
for i = 1:3
    zt = [1 5 10];
    u2 = propTF(U0, L5, lam5, zt(i));
    I2 = abs(u2).^2;
    wz = w0*sqrt(1+(zt(i)/zR)^2);       % 式 (5.27)
    I_an = abs(A0)^2*(w0/wz)^2*exp(-2*(X5.^2+Y5.^2)/wz^2);   % 式 (5.26)
    subplot(1,3,i);
    plot(x5, I_an(M5/2+1,:), '--'); hold on;
    plot(x5, I2(M5/2+1,:), '-');
    grid on; xlabel('x (m)'); ylabel('Irradiance');
    title(sprintf('z = %d m, w(z) = %.2f mm', zt(i), wz*1e3));
    legend('analytic (5.26)', 'propTF', 'FontSize', 8);
    err = max(abs(I2(M5/2+1,:) - I_an(M5/2+1,:)))/max(I_an(M5/2+1,:));
    fprintf('(a) z=%4d m: 剖面最大相对偏差 = %.2e\n', zt(i), err);
end
saveas(gcf, 'out_ex5_5.png');

% ================= 习题 5.6:split-step 等价性(式 5.15) =================
u_single = propTF(u1q, Lq, lamq, 20000);      % 单次 TF 20000 m
u_split = u1q;
for i = 1:10                                  % 10 次 TF,每步 2000 m(每步临界)
    u_split = propTF(u_split, Lq, lamq, 2000);
end
dI = max(max(abs(abs(u_single).^2 - abs(u_split).^2)));
fprintf('\n========== 习题 5.6(split-step 等价性,总距离 20000 m) ==========\n');
fprintf('max |I_single - I_split| = %.3e(浮点误差级 -> 式 (5.15) 成立)\n', dI);
fprintf('结论:拆步不改变采样 regime(总距离才决定);split-step 的用途是\n');
fprintf('在多相位屏(大气湍流)等场景中途插入元件,不是救欠采样。\n');
