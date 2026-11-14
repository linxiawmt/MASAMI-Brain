% ===== 书 2.2 节：有效带宽（式 2.9-2.11 + 表 2.1） =====
% 有效带宽：包含 98% 总频谱功率的（径向）频谱宽度：
%   P(B)/P_T = 0.98,  P(B) = ∫∫_{ρ<B} |F|^2 ρ dρ dθ          （式 2.10）
% 总功率 P_T 用 Parseval 在空域算（式 2.9），数值累积除以解析 P_T，
% 避免 98% 附近的平坦曲线被截断归一化带偏。
% 书表 2.1：rect(x/2w)rect(y/2w) → 5/w；circ(sqrt(x^2+y^2)/w) → 5/w；
%           exp[-π(x^2+y^2)/w^2] → 0.79/w。
% 圆函数半径约定（承接 ch1 勘误 1）：表 2.1 的 circ 行对应"半径 = w"的圆盘
% （Goodman/表 1.3 约定），B·w ≈ 4.99；若按表 1.2 的"半径 = w/2"则是 ≈ 9.98。
% MATLAB/Python 差异：besselj 内置（Python 用 scipy.special.j1）；
% 索引 1-based；cumtrapz/interp1 为基础函数（Python 需自己拼）。
clear; clc; close all;
w   = 1.0;                      % 特征半宽（结果按 B·w 报告）
NR  = 8001;                     % ρ 网格密度
rho = linspace(0, 60/w, NR);

% ---------- rect：|F|^2 非圆对称，做完整 θ 积分（式 2.10） ----------
% 正方形对称性：θ∈[0,π/4] 积分后乘 8
NT = 4096;
th = linspace(0, pi/4, NT);
F2rect = @(fX, fY) (4*w^2 .* sincn(2*w*fX) .* sincn(2*w*fY)).^2;
ring = zeros(size(rho));
for s = 1:128:NR
    e = min(s+127, NR);
    r  = rho(s:e).';                 % 列向量
    fX = r .* cos(th);               % 行 = θ
    fY = r .* sin(th);
    v  = F2rect(fX, fY);
    ring(s:e) = sum((v(:,1:end-1) + v(:,2:end)) / 2, 2) * (th(2)-th(1));
end
ring = ring * 8;
P_rect = [0, cumtrapz(ring .* rho)];        % P(B)

% ---------- circ / gauss：圆对称，径向 1D 即可 ----------
P_circ  = [0, cumtrapz(2*pi .* (jinc(rho, w)).^2        .* rho)];  % 半径 w 圆盘
P_gauss = [0, cumtrapz(2*pi .* (w^2*exp(-pi*w^2*rho.^2)).^2 .* rho)];

P_T = struct('rect', 4*w^2, 'circW', pi*w^2, 'gauss', w^2/2);
B_rect  = interp1(P_rect,  rho, 0.98*P_T.rect)  * w;
B_circ  = interp1(P_circ,   rho, 0.98*P_T.circW) * w;
B_gauss = interp1(P_gauss,  rho, 0.98*P_T.gauss) * w;

fprintf('98%% 有效带宽（数值，径向判据）\n');
fprintf('  rect(x/2w)rect(y/2w) : B·w = %.3f   （书表 2.1：5/w）\n', B_rect);
fprintf('  circ 半径 w (Goodman/表1.3 约定): B·w = %.3f（书表 2.1：5/w ✓）\n', B_circ);
fprintf('  circ 半径 w/2 (表 1.2 定义)     : B·w = %.3f\n', 2*B_circ);
fprintf('  exp(-pi r^2/w^2)     : B·w = %.4f（书表 2.1：0.79/w；解析 %.4f）\n', ...
        B_gauss, sqrt(log(50)/(2*pi)));

% ---------- 画图：左两幅 = rect 频谱 98% 范围（复现图 2.3），右 = 累积功率曲线 ----------
fmax = 6/w; n = 601;
fx1 = linspace(-fmax, fmax, n);
[FX, FY] = meshgrid(fx1, fx1);      % 'xy' 行为，行 ↔ y，列 ↔ x
Fmag = 4*w^2 .* sincn(2*w*FX) .* sincn(2*w*FY);

figure('Position', [80 80 1350 420]);
subplot(1,3,1);
imagesc(fx1, fx1, abs(Fmag)); axis xy; colorbar; hold on;
th_ = linspace(0, 2*pi, 200); plot(B_rect/w*cos(th_), B_rect/w*sin(th_), 'c--', 'LineWidth', 1.5);
xlabel('f_X/w'); ylabel('f_Y/w'); title('|F(f_X,f_Y)|, rect, circle = 98% power');
subplot(1,3,2);
imagesc(fx1, fx1, abs(Fmag).^2); axis xy; colorbar; hold on;
plot(B_rect/w*cos(th_), B_rect/w*sin(th_), 'c--', 'LineWidth', 1.5);
xlabel('f_X/w'); ylabel('f_Y/w'); title('|F|^2 (power spectrum), rect');
subplot(1,3,3);
plot(rho*w, P_rect/P_T.rect, 'LineWidth', 1.2); hold on;
plot(rho*w, P_circ/P_T.circW, 'LineWidth', 1.2);
plot(rho*w, P_gauss/P_T.gauss, 'LineWidth', 1.2);
yline(0.98, 'k:'); plot(B_rect, 0.98, 'o', B_circ, 0.98, 'o', B_gauss, 0.98, 'o', 'MarkerSize', 5);
xlim([0 12]); ylim([0.9 1.001]); grid on;
xlabel('B·w'); ylabel('cumulative power fraction');
title('98% effective bandwidth (Table 2.1)');
legend('rect(x/2w)rect(y/2w)', 'circ, disk R=w (Table 2.1)', 'exp(-\pir^2/w^2)', 'FontSize', 8);
saveas(gcf, 'out_eff_bw.png');
fprintf('已保存 out_eff_bw.png\n');

% ---------- 局部函数（定义同书表 1.2 / ch1 勘误 1 的可靠形式） ----------
function r = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)，u=0 处极限 1
    r = ones(size(u));
    nz = u ~= 0;
    r(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end

function A = jinc(rho, R)
    % 半径 R 圆盘的解析频谱：R*J1(2*pi*R*rho)/rho，rho=0 极限 = 面积
    A = zeros(size(rho));
    nz = rho ~= 0;
    A(nz)  = R .* besselj(1, 2*pi*R*rho(nz)) ./ rho(nz);
    A(~nz) = pi * R^2;
end
