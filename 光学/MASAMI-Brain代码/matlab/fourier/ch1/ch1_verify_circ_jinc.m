% ===== 书 1.5 节：circ 的二维傅里叶变换对验证（原书图 1.3 的 "jinc" 曲面） =====
% 表 1.2 定义 circ(sqrt(x^2+y^2)) 在半径 1/2 内取 1；
% 表 1.3 却给出 circ(r/a) <-> a^2*J1(2*pi*a*rho)/(a*rho)，峰值 pi*a^2（a=1 时为
% pi，正是图 1.3 标注的峰值）。半径 a/2 的圆盘面积是 pi*a^2/4，而频谱峰值应等于
% 面积（中心纵坐标定理），所以表 1.3 这一对应的是"半径 a 的圆盘"（Goodman 约定）。
% 下面用二维 DFT 同时验证两种半径，看解析式到底匹配哪一个。
clear; clc; close all;
N  = 512;           % 每边采样点数
L  = 4.0;           % 空间窗口宽度 [m]
dx = L / N;
n  = (0:N-1) - N/2;
x  = n * dx;
f  = n / (N * dx);  % 频率轴 [cycles/m]
[X, Y] = meshgrid(x, x);   % 圆对称函数下与 numpy 'xy'/'ij' 索引结果一致

Rlist = [1.0 0.5];
results = cell(size(Rlist));       % {f, row, A, err}
for k = 1:numel(Rlist)
    R = Rlist(k);
    C = double(sqrt(X.^2 + Y.^2) < R);          % 半径 R 圆盘（边界取值不影响积分）
    G = fftshift(fft2(ifftshift(C))) * dx^2;    % 二维连续 FT 的 DFT 近似
    row = real(G(N/2+1, :));                    % fY=0 的径向剖面（DC 在 N/2+1）
    A   = jinc(f, R);
    err = max(abs(row - A));
    area = sum(C(:)) * dx^2;                    % 中心纵坐标定理：G(0,0) = 面积
    results{k} = struct('f', f, 'row', row, 'A', A, 'err', err, 'R', R);
    fprintf(['R = %4.2f: 频谱峰值(数值) = %.6f, 解析极限 = %.6f, ' ...
             '圆盘面积(数值) = %.6f, 径向剖面最大误差 = %.3e\n'], ...
            R, row(N/2+1), pi*R^2, area, err);
end

fprintf('\n结论：解析式 R*J1(2*pi*R*rho)/rho 精确匹配"半径 R 的圆盘"。\n');
fprintf('表 1.3（a=1，峰值 pi）对应半径 1 的圆盘；若坚持表 1.2 的半径 1/2 定义，\n');
fprintf('变换应为 0.5*J1(pi*rho)/rho（峰值 pi/4）——两处约定不一致，使用时务必小心。\n');

% ---------- 画径向剖面对比图 ----------
figure('Position', [100 100 1100 400]);
for k = 1:2
    subplot(1, 2, k);
    S = results{k};
    m = abs(S.f) <= 3;
    plot(S.f(m), S.A(m), 'b-', 'LineWidth', 1.5); hold on;
    plot(S.f(m), S.row(m), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'none');
    yl = ylim; plot(xlim, [0 0], 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    ylim(yl); grid on; xlabel('\rho [cycles/m]');
    if S.R == 1
        title(sprintf('R = 1 (peak \\pi, matches Table 1.3 / Fig. 1.3)\nmax error = %.1e', S.err));
    else
        title(sprintf('R = 1/2 (peak \\pi/4, matches Table 1.2 definition)\nmax error = %.1e', S.err));
    end
    legend('analytic R J_1(2\pi \rho)/\rho', '2D DFT numeric', 'Location', 'northwest');
end
sgtitle('Circle function transform (jinc) — radial cut fY = 0');
saveas(gcf, 'out_circ_jinc.png');
fprintf('已保存 out_circ_jinc.png\n');

% ---------- 局部函数 ----------
function A = jinc(rho, R)
    % 半径 R 圆盘的解析频谱：R*J1(2*pi*R*rho)/rho；rho=0 处极限 = 圆盘面积
    A = zeros(size(rho));
    nz = rho ~= 0;
    A(nz)  = R * besselj(1, 2*pi*R*rho(nz)) ./ rho(nz);   % besselj: 第一类贝塞尔
    A(~nz) = pi * R^2;
end
