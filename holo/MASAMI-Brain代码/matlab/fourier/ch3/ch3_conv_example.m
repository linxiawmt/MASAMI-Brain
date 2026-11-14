% ===== 书 3.6 节 + 习题 3.2：两个高斯函数的 FFT 卷积与解析对比 =====
% 参数（书 3.6 节）：wa = 0.3 m、wb = 0.2 m——注释原文 "[exp(-pi) radius]"，
% 意思是 fa(wa) = exp(-pi) ≈ 0.043 的半径（不是 1/e 半径！）；L = 2、M = 200。
% 支撑判据（式 2.27）：Da ≈ 0.9、Db ≈ 0.6，和 1.5 < L = 2 -> 卷绕可忽略。
% 书的技巧：fa、fb 都居中排列，FFT 前【不移位】——两个函数携带同一线性
% 相位，乘积相位对应"整周期平移"，模 L 不可见；但逆变换结果停在"移位
% 排列"（峰值在 index 1），所以要 fftshift 居中。只移位其中一个函数，
% 卷积结果就会整体错位 L/2（本脚本有反例演示）。
% MATLAB/Python 差异：.* 与 .^ 的点号在 numpy 里不需要——Python 的 * 与
% ** 天然逐元素，矩阵乘法反而要另写 @ / np.matmul，方向恰好相反。
clear; clc; close all;

wa = 0.3;                   % 高斯 a 宽度 [exp(-pi) 半径] (m)
wb = 0.2;                   % 高斯 b 宽度 [exp(-pi) 半径] (m)
L  = 2;                     % 边长 (m)
M  = 200;                   % 样本数
dx = L/M;                   % 采样间隔 (m)

x  = -L/2:dx:L/2-dx;        % x 坐标（式 2.21）
fa = exp(-pi*(x.^2)/wa^2);  % 高斯 a（注意 .^ 是逐元素平方）
fb = exp(-pi*(x.^2)/wb^2);  % 高斯 b

% ---------- FFT 卷积（书 3.6 节，输入不移位） ----------
Fa = fft(fa);               % 变换 fa
Fb = fft(fb);               % 变换 fb
F0 = Fa.*Fb;                % 逐点相乘（.* 不是 *，忘点号是最常见错误！）
f0 = ifft(F0)*dx;           % 逆变换；卷积定标只乘一次 dx（不是 dx^2）
f  = fftshift(f0);          % 居中（书图 3.9(b)）

% ---------- 解析卷积（习题 3.2） ----------
% exp(-pi*x^2/w^2) <-> w*exp(-pi*w^2*fx^2)（表 1.3 高斯对）。频域相乘：
% wa*wb*exp(-pi*(wa^2+wb^2)*fx^2)，逆变换回去：
%   c(x) = wa*wb/sqrt(wa^2+wb^2) * exp(-pi*x^2/(wa^2+wb^2))
c_an = wa*wb/sqrt(wa^2+wb^2) * exp(-pi*x.^2/(wa^2+wb^2));

fprintf('解析卷积峰值 = %.6f；FFT 卷积峰值 = %.6f\n', max(c_an), max(f));
fprintf('FFT 卷积与解析最大差异 = %.2e（高斯衰减极快：离散求和与连续积分之差达机器精度）\n', ...
        max(abs(f - c_an)));

% ---------- 反例：只移位一个函数 -> 卷积整体错位 ----------
fa_s    = ifftshift(fa);                        % 只把 fa 挪到移位排列
f_wrong = fftshift(real(ifft(fft(fa_s).*fft(fb))*dx));
[~, iw] = max(f_wrong);
fprintf('只移位 fa：卷积峰出现在 x = %.2f m（应为 0 -> 整体错位 L/2）\n', x(iw));

% ---------- 画图：复现书图 3.9 + 解析对比 ----------
figure('Position', [80 80 950 380]);
subplot(1,2,1);
plot(x, fa, x, fb, '--'); grid on;
xlabel('x (m)'); title('functions (Fig. 3.9a)');
legend('f_a (w_a=0.3)', 'f_b (w_b=0.2)', 'FontSize', 8);
subplot(1,2,2);
plot(x, f, x, c_an, ':'); grid on;
xlabel('x (m)'); title('convolution: FFT vs analytic (Fig. 3.9b + Ex. 3.2)');
legend('FFT', 'analytic', 'FontSize', 8);
saveas(gcf, 'out_conv_example.png');
fprintf('已保存 out_conv_example.png\n');
