function [phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0)
% ===== 书 Schmidt PM199 Listing 9.3 =====
% FT 法 + 次谐波（subharmonics）生成大气湍流相位屏（Lane et al. 方法）：
%   phz = phz_hi + phz_lo；
%   phz_hi：Listing 9.2 的普通 FT 屏（高频部分）；
%   phz_lo：低频补偿屏（式 9.81）——对 p = 1,2,3 三层"次谐波网格"
%   （间距 1/(3^p * L)，L = N*delta）各取 3×3 = 9 个最低频分量，
%   按 PSD 平方根抽复高斯系数后直接求和（频率太少，不必用 FFT）。
% 低频网格把 FT 屏欠采样的倾斜/低阶模式补回来，使大间隔结构函数
% 与理论 6.88*(r/r0)^(5/3) 吻合（图 9.3）。
%
% 用法：[phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0)
% 注意：内部用 randn，未播种；需要可复现时在调用前 rng(seed)。

D = N*delta;                             % 屏的边长 [m]
% high-frequency screen from FFT method
phz_hi = ft_phase_screen(r0, N, delta, L0, l0);
% spatial grid [m]
[x y] = meshgrid((-N/2 : N/2-1) * delta);
% initialize low-freq screen
phz_lo = zeros(size(phz_hi));

% loop over frequency grids with spacing 1/(3^p*L)
for p = 1:3
    % setup the PSD
    del_f = 1 / (3^p*D);                 % frequency grid spacing [1/m]
    fx = (-1 : 1) * del_f;               % frequency grid [1/m]（3×3 个频点）
    [fx fy] = meshgrid(fx);
    [th f] = cart2pol(fx, fy);           % polar grid
    fm = 5.92/l0/(2*pi);                 % inner scale frequency [1/m]
    f0 = 1/L0;                           % outer scale frequency [1/m]

    % modified von Karman atmospheric phase PSD
    PSD_phi = 0.023*r0^(-5/3) * exp(-(f/fm).^2) ...
        ./ (f.^2 + f0^2).^(11/6);
    PSD_phi(2,2) = 0;                    % 网格中心 (0,0) 频点置零

    % random draws of Fourier coefficients（式 9.79 同款缩放）
    cn = (randn(3) + i*randn(3)) ...
        .* sqrt(PSD_phi)*del_f;
    SH = zeros(N);
    % loop over frequencies on this grid（9 个频点逐个叠加平面波）
    for ii = 1:9
        SH = SH + cn(ii) ...
            * exp(i*2*pi*(fx(ii)*x+fy(ii)*y));
    end
    phz_lo = phz_lo + SH;                % accumulate subharmonics
end
phz_lo = real(phz_lo) - mean(real(phz_lo(:)));   % 取实部并去掉均值（piston）
end
