%% ===== 书 Schmidt PM199 第 1 章 =====
%  几何光学 vs 波动光学：不透明圆盘（水滴模型）挡光后的光强分布。
%  几何光学给出"绝对阴影"；波动光学（第 1 章的菲涅耳衍射积分，式 1.57）给出
%  阴影边缘的衍射条纹，以及阴影中心残留的亮斑（泊松斑）。
%  重写自用户旧稿 几何光学和波动光学对比.py：旧版用 1-(D/λz)·sinc 的启发式公式示意，
%  本版改为按式 (1.57) 的传递函数法真实传播，物理量纲正确。
clear; clc; close all;

%% 1. 参数
N   = 1024;
dx  = 2e-6;           % 2 um 像素
wvl = 532e-9;
k   = 2*pi/wvl;
z   = 5e-3;           % 传播 5 mm（传递函数法窗口内：z < 4 dx^2/λ ≈ 30 mm）

x = ((1:N) - N/2 - 0.5) * dx;
[X, Y] = meshgrid(x, x);
R = sqrt(X.^2 + Y.^2);

D = 100e-6;           % 圆盘（水滴）直径 100 um

%% 2. 两个模型
% 几何光学：圆盘内光强为 0，其余为 1（绝对阴影）
I_geom = double(R > D/2);

% 波动光学：U1 = 1 - circ（平面波减圆盘），用式(1.57)的 TF 形式传播
U1 = double(R > D/2);
fx = ((1:N) - N/2 - 1) / (N*dx);
[FX, FY] = meshgrid(fx, fx);
H = exp(1i*k*z) * exp(-1i*pi*wvl*z*(FX.^2 + FY.^2));
U_ft = fftshift(fft2(fftshift(U1)));
U2 = ifftshift(ifft2(ifftshift(U_ft .* H)));
I_wave = abs(U2).^2;

%% 3. 数值验证
fprintf('===== 几何 vs 波动数值验证 =====\n');
ctr = N/2 + 1;
% (1) 泊松斑：几何光学预言盘心 I=0；波动光学预言轴上振幅恢复到入射值（I=1）。
%     网格离散下泊松斑尺度 sqrt(lambda*z)=1.6 um 小于像素，中心值是它的像素平均残影。
fprintf('(1) 盘心光强：几何 = 0，波动数值 = %.3f（理想连续值 1）\n', I_wave(ctr, ctr));
% (2) 阴影边缘：几何光学的跳变 vs 波动光学的条纹（中心行剖面）
row_g = I_geom(ctr, :);  row_w = I_wave(ctr, :);
out_sh = abs(x) > D/2 & abs(x) < D/2 + 50e-6;   % 阴影外 50 um 环带
fprintf('(2) 阴影外 50 um 环带内 波动光强范围 = [%.3f, %.3f]（几何恒为 1，衍射产生明暗条纹）\n', ...
    min(row_w(out_sh)), max(row_w(out_sh)));
% (3) 能量守恒：TF 传播前后总能量（|U|^2 求和）应基本不变
fprintf('(3) 传播前 sum|U1|^2 = %.1f，传播后 sum|U2|^2 = %.1f\n', ...
    sum(U1(:).^2), sum(I_wave(:)));

%% 4. 绘图
figure('Name','几何光学 vs 波动光学','Position',[80 80 1200 460]);

subplot(1,2,1);
imagesc(x*1e6, x*1e6, I_wave); axis xy; axis image; colormap inferno; colorbar;
caxis([0 1]);
% 用参数方程画虚线圆周（不用 Image Processing Toolbox 的 viscircles）
th = linspace(0, 2*pi, 181);
plot(D/2*1e6*cos(th), D/2*1e6*sin(th), 'w--', 'LineWidth', 0.8);
xlim([-200 200]); ylim([-200 200]);
title(sprintf('波动光学：圆盘后 z=%d mm 的光强分布', z*1e3));
xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(1,2,2);
plot(x*1e6, row_g, 'r--', 'LineWidth', 1.5); hold on;
plot(x*1e6, row_w, 'b-', 'LineWidth', 1.5);
xline(-D/2*1e6, 'Color', [0.5 0.5 0.5]); xline(D/2*1e6, 'Color', [0.5 0.5 0.5]);
xlim([-150 150]); ylim([-0.05 1.15]); grid on;
xlabel('x (\mum)'); ylabel('相对光强');
legend('几何光学（绝对阴影）', '波动光学（菲涅耳衍射）', 'Location', 'southeast');
title('中心行剖面：两种模型的对比');
