% ===== 书 Schmidt PM199 Listing 9.5 =====
% 例题（9.5.1 节）：点源经 50 km、Cn2 = 1e-16 m^-2/3 的湍流路径传到
% 0.5 m 口径望远镜 —— 建立几何、算大气参数、用约束优化解 11 块相位屏
% 各自的 r0（分层矩匹配，式 9.65-9.75）。
% 关键结果（原书正文）：r0,sw = 12.7 cm、sigma_chi^2,sw = 0.436。
% 依赖：Optimization Toolbox 的 fmincon（带边界约束的最小化）；
%       Python 版用 scipy.optimize.lsq_linear（本例是无约束rank的
%       有界线性最小二乘，两者等价，见 python/waveprop/ch9）。
% 运行前提：ft2.m / ift2.m 在 path 上（本文件本身不用，但同目录脚本需要）。

% determine geometry
D2 = 0.5;            % diameter of the observation aperture [m]
wvl = 1e-6;          % optical wavelength [m]
k = 2*pi / wvl;      % optical wavenumber [rad/m]
Dz = 50e3;           % propagation distance [m]

% use sinc to model pt source（点源的 sinc-Gaussian 模型，见 Listing 9.7）
DROI = 4 * D2;           % diam of obs-plane region of interest [m]
D1 = wvl*Dz / DROI;      % width of central lobe [m]（sinc 中央瓣全宽 2*D1）
R = Dz;                  % wavefront radius of curvature [m]

% atmospheric properties
Cn2 = 1e-16;         % structure parameter [m^-2/3], constant（沿路径常数）

% SW and PW coherence diameters [m]（式 9.43 / 9.42，常数 Cn2 的解析积分）
r0sw = (0.423 * k^2 * Cn2 * 3/8 * Dz)^(-3/5);   % 球面波（点源）
r0pw = (0.423 * k^2 * Cn2 * Dz)^(-3/5);         % 平面波

p = linspace(0, Dz, 1e3);
% log-amplitude variance（式 9.64 的数值积分，中点求和）
rytov = 0.563 * k^(7/6) * sum(Cn2 * (1-p/Dz).^(5/6) ...
    .* p.^(5/6) * (p(2)-p(1)));

% screen properties
nscr = 11;                          % number of screens（11 个平面 / 10 次部分传播）
A = zeros(2, nscr);                 % matrix（式 9.75 的 2×nscr 矩阵）
alpha = (0:nscr-1) / (nscr-1);      % 各屏的相对位置 z_i / Dz（均匀）
A(1,:) = alpha.^(5/3);              % 第一行：r0,sw 方程的权重 alpha^(5/3)
A(2,:) = (1 - alpha).^(5/6) .* alpha.^(5/6);   % 第二行：Rytov 方程权重
b = [r0sw.^(-5/3); rytov/1.33*(k/Dz)^(5/6)];   % 式 9.72 / 9.74 移项后的左端

% initial guess（每屏 r0 初值 = nscr/3 * r0sw，取其 -5/3 次幂作未知量）
x0 = (nscr/3*r0sw * ones(nscr, 1)).^(-5/3);

% objective function（min ||A*X - b||^2，X = [r0_i^(-5/3)]）
fun = @(X) sum((A*X(:) - b).^2);

% constraints：下界 0（r0 -> 无穷 = 无湍流）；上界使每屏对 Rytov 数的
% 贡献 <= rmax = 0.1（Martin & Flatté 指南：单屏/单步 Rytov 数要小）
x1 = zeros(nscr, 1);
rmax = 0.1;                         % maximum Rytov number per partial prop
x2 = rmax/1.33*(k/Dz)^(5/6) ./ A(2,:);
x2(A(2,:)==0) = 50^(-5/3)           % 权重为 0 的屏（路径两端）给个宽松上界（r0 = 50 m）

[X,fval,exitflag,output] ...
    = fmincon(fun,x0,[],[],[],[],x1,x2)

% check screen r0s
r0scrn = X.^(-3/5)
r0scrn(isinf(r0scrn)) = 1e6;        % X = 0（下界）的屏 -> r0 = 无穷，记 1e6 m

% check resulting r0sw & rytov（回代 vs 目标，应逐项吻合）
bp = A*X(:); [bp(1)^(-3/5) bp(2)*1.33*(Dz/k)^(5/6)]
[r0sw rytov]
