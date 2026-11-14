% ===== 书 4.4.1/4.4.2 节：RS 传递函数 (4.19) vs 菲涅耳传递函数 (4.24) =====
% 对比两张传递函数：菲涅耳近似的可信频率范围 + |f|>1/lambda 的倏逝波截止。
% 相位差推导：
%   arg(H_RS/H_F) = kz(sqrt(1-(lam f)^2) - 1) + pi*lam*z*f^2
%     ≈ -pi*lam^3*z*f^4/4 （式 (4.19) 根号展开到 (lam f)^4 阶时，二阶项与
%       式 (4.24) 的 -pi*lam*z*f^2 正好抵消——"二阶精度"近似的含义）
% MATLAB/Python 差异：MATLAB 的 sqrt 对负实数自动返回复数，Python 的
%   np.sqrt 给 nan、需输入复数 (a+0j)；画图均只取幅度/相位。
clear; clc; close all;

lam = 0.633e-6;          % 波长 (m)，He-Ne
z   = 50;                % 传播距离 (m)（书 4.5 节例子的距离）
k   = 2*pi/lam;

f    = linspace(0, 1.5/lam, 4000);   % 频率轴：以 1/lambda 为单位标定
prop = f < 1/lam;                    % 传播波分量（式 4.19 下方的条件）

% ---------- |H_RS|：传播区幅度 1，倏逝区指数衰减 ----------
mag_RS = zeros(size(f));
mag_RS(prop)  = 1;                                   % |exp(j 相位)| = 1
mag_RS(~prop) = exp(-k*z*sqrt((lam*f(~prop)).^2 - 1));   % 倏逝衰减

% ---------- 相位差（只定义在传播区） ----------
fp    = f(prop);
dphi  = abs(k*z*(sqrt(1-(lam*fp).^2) - 1) + pi*lam*z*fp.^2);
fa    = fp(fp > 1e4);                        % 渐近线只在大 f 画
asym  = pi*lam^3*z*fa.^4/4;
f_01  = (0.4/(pi*lam^3*z))^0.25;             % |dphi| = 0.1 rad 的频率

fprintf('1/lambda（传播波截止）= %.3e cyc/m\n', 1/lam);
fprintf('|dphi|=0.1 rad 的频率 = %.3e cyc/m = %.4f /lambda\n', f_01, f_01*lam);
fprintf(['例：f = 0.01/lambda 处相位差 = %.3f rad' ...
         '（f^4 增长：f 每大 10 倍，误差大 1 万倍）\n'], ...
         pi*lam^3*z*(0.01/lam)^4/4);

figure('Position', [80 80 1050 420]);
subplot(1,2,1);
h1 = semilogy(f*lam, mag_RS); hold on;
yline(1, '--', 'Color', 'crimson');
xline(1, ':', 'Color', [0.5 0.5 0.5]);
ylim([1e-14 3]); xlabel('f_X \lambda'); ylabel('magnitude');
title('RS transfer function: evanescent cut at |f| = 1/\lambda');
legend([h1], {'|H_{RS}|  Eq. (4.19)  (|H_F| = 1 flat)'}, 'FontSize', 8);
grid on;

subplot(1,2,2);
h2 = loglog(fp*lam, dphi); hold on;
h3 = loglog(fa*lam, asym, ':');
yline(0.1, '--', 'Color', 'crimson');
xline(f_01*lam, ':', 'Color', [0.5 0.5 0.5]);
xlabel('f_X \lambda'); ylabel('phase difference (rad)');
title('Fresnel vs RS phase error (z = 50 m)');
legend([h2 h3], {'|arg(H_{RS}/H_F)|', '\pi\lambda^3 z f_X^4/4'}, 'FontSize', 8);
grid on;
saveas(gcf, 'out_transfer_compare.png');
fprintf('已保存 out_transfer_compare.png\n');
