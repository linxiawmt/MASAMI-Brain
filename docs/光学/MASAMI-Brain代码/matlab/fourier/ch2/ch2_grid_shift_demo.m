% ===== 书 2.4 节：坐标、索引、居中与移位（式 2.21/2.22，复现图 2.4 / 2.5） =====
% 三种排列（书图 2.4）：
%   analytic ：指标带正负（m = -M/2..M/2-1），纸面公式用；
%   centered ：零坐标在 index M/2+1（MATLAB 1-based）；
%   shifted  ：index 1 = 零坐标，FFT 前必须转成它。
% centered <-> shifted 用 fftshift / ifftshift（偶数 N 两者相同，奇数 N 不同！）。
% MATLAB/Python 对照：居中零坐标 MATLAB (N/2+1, M/2+1) <-> Python (N//2, M//2)；
% MATLAB (i,j) = 行,列：j ↔ x 坐标、i ↔ y 坐标（与笛卡尔 (x,y) 顺序相反）。
clear; clc; close all;

% ---------- 1D：N=8 小网格，rect 半宽 1.5（3 个样本为 1） ----------
N = 8; L = 8.0; dx = L/N;
n = (0:N-1) - N/2;
x = n * dx;                       % 式 (2.21)：-L/2 : Δx : L/2-Δx
f = n / (N*dx);                   % 式 (2.22)：-1/(2Δx) : 1/L : 1/(2Δx)-1/L
g_centered = rectn(x/3);          % 居中排列：|x|<1.5 → x=-1,0,1
g_shifted  = ifftshift(g_centered);

fprintf('N=8, L=8, dx=1 的排列（对应书图 2.4）\n');
fprintf('MATLAB 1-based index : %s\n', mat2str(1:N));
fprintf('coords x (式 2.21)   : %s\n', mat2str(x));
fprintf('g 居中 centered      : %s\n', mat2str(g_centered));
fprintf('g 移位 shifted(FFT用): %s\n', mat2str(g_shifted));
fprintf('freq 坐标 f (式 2.22): %s\n', mat2str(f, 4));
fprintf('roundtrip fftshift(ifftshift(g)) == g : %d\n', ...
        isequal(fftshift(g_shifted), g_centered));
fprintf(['注意：偶数 N 时 fftshift 与 ifftshift 等价；奇数 N 不等价' ...
         '（见笔记坑与陷阱）\n\n']);

% ---------- 2D：象限交换（书图 2.5） ----------
N2 = 8;
[X, Y] = meshgrid(x, x);          % 行 ↔ y，列 ↔ x（与 Python 'xy' 一致）
A_centered = double(sqrt(X.^2 + Y.^2) < 1.5);   % 居中圆盘
A_shifted  = ifftshift(A_centered);
[iz, jz] = find(A_centered > 0 & abs(X) < 1e-12 & abs(Y) < 1e-12);
fprintf('2D 居中排列零坐标位置 (行,列) = (%d, %d)；即 (N/2+1, M/2+1)\n', iz, jz);
fprintf('2D 移位排列零坐标位置 (行,列) = (1, 1)\n\n');

% ---------- 画图：上排复现图 2.4 三种 1D 排列，下排复现图 2.5 象限交换 ----------
figure('Position', [60 60 1150 760]);

subplot(2,3,1);
plot(x, g_centered, 'o-', 'MarkerSize', 5, 'Color', 'crimson'); grid on;
xlim([-4.5 3.5]); ylim([-0.2 1.25]);
title('centered: zero coord at index N/2+1');
xlabel('x');

subplot(2,3,2);
plot(0:N-1, g_shifted, 's-', 'MarkerSize', 5, 'Color', [0 .55 .55]); grid on;
ylim([-0.2 1.25]);
title('shifted for FFT: zero coord at index 1');
xlabel('array index (1-based)');

subplot(2,3,3);
G = fft(g_shifted) * dx;
plot(f, abs(fftshift(G)), 'o-', 'MarkerSize', 5, 'Color', 'slateblue'); grid on;
title('|FFT \times \Deltax| back on centered freq axis');
xlabel('f_X');

subplot(2,3,4);
imagesc(qmap(N2)); axis image off;   % 居中：3,4 / 2,1
title('2D centered (quadrant No.), [0,0] at (N/2+1, M/2+1)');
subplot(2,3,5);
imagesc(ifftshift(qmap(N2))); axis image off;  % 移位：1,2 / 4,3
title('2D shifted for FFT (quadrant No.), [0,0] at (1,1)');

subplot(2,3,6);
plot(abs(G), 's-', 'MarkerSize', 5, 'Color', [1 .55 0]); grid on;
xline(N2/2 + 0.5, ':');
title('why fftshift exists: raw fft ordering');
xlabel('array index (1-based)');

saveas(gcf, 'out_arrangements.png');
fprintf('已保存 out_arrangements.png\n');

% ---------- 局部函数 ----------
function r = rectn(t)
    r = zeros(size(t));
    r(abs(t) < 0.5) = 1;
    r(abs(t) == 0.5) = 0.5;
end

function q = qmap(N)
    % 象限编号（书图 2.5 的居中排列：左上 3、右上 4、左下 2、右下 1）
    q = zeros(N);
    q(1:N/2, 1:N/2)         = 3;
    q(1:N/2, N/2+1:end)     = 4;
    q(N/2+1:end, 1:N/2)     = 2;
    q(N/2+1:end, N/2+1:end) = 1;
end
