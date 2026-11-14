function phz = ft_phase_screen(r0, N, delta, L0, l0)
% ===== 书 Schmidt PM199 Listing 9.2 =====
% FT 法生成大气湍流相位屏（McGlamery 方法，式 9.76-9.80）：
%   1) 在频域按修正 von Karman 相位 PSD（式 9.51 的 f 换算版，式 9.52）
%      的平方根 * 生成复高斯随机系数 cn（式 9.79：E{|cn|^2} = PSD*df^2）；
%   2) 逆 FT 合成相位屏；real(IFT) 的实部/虚部是两块不相关的屏，
%      这里取实部、丢弃虚部。
% 缺陷（书 9.3 节，图 9.3）：频域网格最低频 = 1/(N*delta)，低阶模式
% （倾斜等）欠采样 -> 大间隔结构函数偏低；补救见 ft_sh_phase_screen。
%
% 用法：phz = ft_phase_screen(r0, N, delta, L0, l0)
%   r0:    Fried 参数 / 相干直径 [m]（本屏的强度）
%   N:     每边网格点数
%   delta: 网格间距 [m]
%   L0:    外尺度 [m]（传 inf 得 Kolmogorov 谱）
%   l0:    内尺度 [m]（传 0 得无内尺度谱）
% 注意：函数内部用 randn，未播种；需要可复现时在调用前 rng(seed)。

% setup the PSD
del_f = 1/(N*delta);                 % frequency grid spacing [1/m]
fx = (-N/2 : N/2-1) * del_f;         % frequency grid [1/m]
[fx fy] = meshgrid(fx);
[th f] = cart2pol(fx, fy);           % polar grid
fm = 5.92/l0/(2*pi);                 % inner scale frequency [1/m]（式 9.18：kappa_m = 5.92/l0，f = kappa/2pi）
f0 = 1/L0;                           % outer scale frequency [1/m]（kappa_0 = 2*pi/L0）

% modified von Karman atmospheric phase PSD（式 9.51 -> 式 9.52：0.0229 ~= 0.023）
% l0 = 0 时 fm = inf、指数 = 1；L0 = inf 时 f0 = 0 —— 退化为 Kolmogorov 谱
PSD_phi = 0.023*r0^(-5/3) * exp(-(f/fm).^2) ...
    ./ (f.^2 + f0^2).^(11/6);
PSD_phi(N/2+1,N/2+1) = 0;            % 零频分量置零（整体 piston 无意义）

% random draws of Fourier coefficients（式 9.79/9.80：E{|cn|^2} = PSD*df^2）
% (randn + i*randn) 的模方期望是 2 倍，取 real(IFT) 又带来 1/2，恰好相消
cn = (randn(N) + i*randn(N)) .* sqrt(PSD_phi)*del_f;

% synthesize the phase screen（IFT 用 delta_f = 1 的"原始傅里叶级数"合成）
phz = real(ift2(cn, 1));
end
