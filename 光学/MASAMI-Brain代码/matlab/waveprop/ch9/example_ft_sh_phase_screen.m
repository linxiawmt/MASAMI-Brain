% ===== 书 Schmidt PM199 Listing 9.4 =====
% ft_sh_phase_screen 的使用示例：生成一块 r0 = 10 cm 的大气相位屏
% （原书图 9.4 即本例的一块典型屏幕）。
% 屏幕参数：边长 2 m、N = 256、外尺度 100 m、内尺度 1 cm。
% 运行前提：本文件与 ft_sh_phase_screen.m / ft_phase_screen.m /
% ift2.m / ft2.m 在同一目录（或均在 MATLAB path 上）。

D = 2;       % length of one side of square phase screen [m]
r0 = 0.1;    % coherence diameter [m]（Fried 参数）
N = 256;     % number of grid points per side
L0 = 100;    % outer scale [m]
l0 = 0.01;   % inner scale [m]

delta = D/N; % grid spacing [m]

% spatial grid
x = (-N/2 : N/2-1) * delta;
y = x;

% generate a random draw of an atmospheric phase screen
[phz_lo phz_hi] ...
    = ft_sh_phase_screen(r0, N, delta, L0, l0);
phz = phz_lo + phz_hi;

% 可视化（对照原书图 9.4）
figure;
imagesc(x, y, phz); axis image; colorbar;
xlabel('x [m]'); ylabel('y [m]'); title('大气相位屏（FT + 次谐波法）');
