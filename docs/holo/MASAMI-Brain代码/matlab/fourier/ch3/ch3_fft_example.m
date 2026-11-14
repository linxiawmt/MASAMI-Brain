% ===== 书 3.2-3.5 节：fft_example 脚本（建向量 -> 移位 -> FFT 定标 -> 解析对比） =====
% 书把这几节写成一个逐步加长的脚本 fft_example.m（原书第 1-58 行），本文件
% 按同样的参数与步骤合并复现书图 3.2-3.8，并附数值检查。
% 参数（书 3.2 节）：矩形半宽 w = 0.055 m（全宽 0.11 m）、边长 L = 2 m、
% 样本数 M = 200 -> 采样间隔 dx = 0.01 m。
% MATLAB/Python 差异：
%   - MATLAB 索引 1-based：居中排列的零坐标在 index M/2+1（本例 101）；
%     Python 0-based 在 index M//2（本例 100）；
%   - 书用 fftshift 做 FFT 前移位：偶数 M 时与 ifftshift 等价；奇数 M 必须
%     用 ifftshift 配对（书 3.3 节原话也这么提醒）；
%   - plot 的 '-o'、axis([xmin xmax ymin ymax]) 对应 Python 的
%     plt.plot(..., 'o-') 与 plt.xlim/plt.ylim。
clear; clc; close all;

% ---------- 3.2 创建向量 ----------
w  = 0.055;                 % 矩形半宽 (m)。故意不取 0.05：见 3.1 节"奇数样本"特性
L  = 2;                     % 向量边长 (m)
M  = 200;                   % 样本数
dx = L/M;                   % 采样间隔 (m)
x  = -L/2:dx:L/2-dx;        % 坐标向量（式 2.21）：-1 : 0.01 : 0.99，共 200 个
f  = rect(x/(2*w));         % 采样矩形：|x| <= 0.055 的样本为 1（11 个）
                            % rect.m 是 3.1 节建的函数（同文件夹即可调用）

fprintf('x 范围: [%.2f, %.2f]，共 %d 个样本\n', x(1), x(end), numel(x));
fprintf('f==1 的样本数 = %d（全宽 0.11 m / dx 0.01 m = 11，奇数）\n', sum(f == 1));
fprintf('坐标 x=0 的 index = %d（居中排列 M/2+1 = %d）\n', ...
        find(abs(x) < 1e-12), M/2 + 1);

% ---------- 3.3 为 FFT 移位 ----------
f0 = fftshift(f);           % 居中 -> 移位排列（零坐标样本挪到 index 1）
% 书图 3.5：移位后中心样本在 index 1，5 个样本在"右"（index 2..6），
% 其余 5 个在"左"——挂在向量末尾（index 197..200）。
% 不移位直接 fft 的后果：相当于对平移了 -L/2 的函数做变换，频域多出一线性
% 相位因子 exp(+j*2*pi*f*L/2)（位移定理）——幅度不变、相位全错。

% ---------- 3.4 计算 FFT 并定标 ----------
F0 = fft(f0) * dx;          % FFT 并乘 dx 回连续量纲（式 2.19 + 黎曼和因子）
F  = fftshift(F0);          % 居中排列
fx = -1/(2*dx) : 1/L : 1/(2*dx) - 1/L;    % 频率坐标（式 2.22）：-50 : 0.5 : 49.5

% ---------- 3.5 与解析结果对比 ----------
% rect(x/(2w)) 的解析频谱（表 1.3 第一行，a = 2w）：F(fx) = 2w sinc(2w fx)。
% 书用 MATLAB 内置 sinc（信号处理工具箱函数，定义与表 1.2 一致、等于
% np.sinc）；这里自带局部函数 sincn，免依赖工具箱。
F_an = 2*w*sincn(2*w*fx);

% ---------- 数值检查 ----------
[~, ic] = max(abs(F));      % 峰值（DC）位置
fprintf('DC 样本 |F(fx=0)| = %.6f（= 矩形面积 2w = %.2f，中心纵坐标定理）\n', ...
        abs(F(ic)), 2*w);
in_main = abs(fx) <= 25;
fprintf(['|FFT|-|解析| 最大差异：全频带 %.4f；|fx|<=25 内 %.2e\n'], ...
        max(abs(abs(F) - abs(F_an))), ...
        max(abs(abs(F(in_main)) - abs(F_an(in_main)))));
fprintf('  （边缘差异 = 周期延拓带来的混叠，书 3.5 节原话"FFT 边缘略高"）\n');
fprintf('sinc 零点间隔 1/(2w) 折合样本数：L/(2w) = %.2f（书说"约每 18 个点"）\n', L/(2*w));
fprintf('峰值处相位 = %.2e（实偶函数 -> 相位只有 0 与 ±pi，±pi 同一相位 mod 2pi）\n', ...
        angle(F(ic)));

% ---------- 画图：复现书图 3.3-3.8 ----------
figure('Position', [60 60 1250 680]);
subplot(2,3,1);                              % 书图 3.3（放大看采样矩形）
plot(x, f, '-o', 'MarkerSize', 4); axis([-0.2 0.2 0 1.5]); grid on;
xlabel('x (m)'); title('f = rect(x/(2w)) samples (Fig. 3.3)');
subplot(2,3,2);                              % 书图 3.4
plot(f, '-o', 'MarkerSize', 4); axis([80 120 0 1.5]); grid on;
xlabel('index'); title('vs index: center at 101 = M/2+1 (Fig. 3.4)');
subplot(2,3,3);                              % 书图 3.5
plot(f0, '-o', 'MarkerSize', 4); axis([0 200 0 1.5]); grid on;
xlabel('index'); title('after fftshift: center at index 1 (Fig. 3.5)');
subplot(2,3,4);                              % 书图 3.6(a)
plot(abs(F0)); grid on; xlabel('index');
title('|F_0| vs index (Fig. 3.6a)');
subplot(2,3,5);                              % 书图 3.6(b)
plot(angle(F0)); grid on; xlabel('index');
title('phase vs index: only 0 and ±\pi (Fig. 3.6b)');
subplot(2,3,6);                              % 书图 3.8(a)
plot(fx, abs(F), fx, abs(F_an), ':'); grid on;
xlabel('fx (cyc/m)'); title('centered |F| vs analytic (Fig. 3.8a)');
legend('discrete', 'analytic', 'FontSize', 8);
saveas(gcf, 'out_fft_example.png');
fprintf('已保存 out_fft_example.png\n');

% ---------- 局部函数 ----------
function s = sincn(u)
    % 归一化 sinc：sin(pi*u)/(pi*u)，u=0 处极限 1（表 1.2；np.sinc 同定义）
    s = ones(size(u));
    nz = u ~= 0;
    s(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end
