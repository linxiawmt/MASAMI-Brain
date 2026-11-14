% ===== 书 2.5.1 · Table 2.2 Fresnel_zone_plate.m（规范化重写） =====
% 画正弦菲涅耳波带片（FZP，点源全息图，式 (2.5-3)/(2.5-4)）：
%     FZP(x,y;z0) = A + B * sin( sigma * ((x-x0)^2 + (y-y0)^2) ),  sigma = 1/z0
% 复现原书图 2.19（z0 与 2z0 的深度依赖）与图 2.20（横向位置 (x0,y0) 依赖）。
% 所有尺度沿用原书"任意单位"（sigma=1/z0，x 步长 0.1，范围 ±12.8）；
% 显示 127*(1+imag(fFZP)) 即 127*(1+sin(sigma*r'^2))，灰度 0~254。
%
% 差异提示：
% 1) 原书用 input() 让你现场输入 z0（2~10）与 x0（-8~8），本版固定四组参数
%    一次性复现图 2.19 与 2.20；想改值改 z1 / x0 即可。
% 2) 原书 image(R) 的行方向与数学坐标相反，meshgrid 后注意 (row,col)=(y,x)。

clear
ROWS = 256; COLS = 256;
s = (0:255)*0.1 - 12.8;                 % 原书 x=-12.8:0.1:12.7（步长 0.1，共 256 点）
[X, Y] = meshgrid(s, s);                % X(row,col) 对应 y 方向，与原书 fFZP(r,c) 一致

z1 = 4;                                 % 基准深度（原书 input z0，2~10 内取值）
x0 = 5;                                 % 图 2.20 的横向偏移（原书 y0=-x0）

% ---- 图 2.19 复现：同一点源（on-axis），深度 z0 与 2*z0 ----
sigma1 = 1/z1;
fFZP1 = exp(1j*sigma1*(X.^2 + Y.^2));           % 原书：exp(j*sigma*(x-x0)^2+j*sigma*(y-y0)^2)
R1 = 127*(1 + imag(fFZP1));                     % 归一化后取虚部 -> 1+sin(phase)

sigma2 = 1/(2*z1);
fFZP2 = exp(1j*sigma2*(X.^2 + Y.^2));
R2 = 127*(1 + imag(fFZP2));

% ---- 图 2.20 复现：同深度，横向位置从 (0,0) 挪到 (x0,y0) ----
sigma3 = 1/z1;
fFZP3 = exp(1j*sigma3*((X-x0).^2 + (Y+x0).^2)); % 原书 y0=-x0
R3 = 127*(1 + imag(fFZP3));

% ---- 可视化 ----
colormap(gray(255))
figure(1)
subplot(2,2,1); image(R1); axis square off
title(sprintf('FZP, on-axis, z0 = %g (fringes coarser when deeper)', z1))
subplot(2,2,2); image(R2); axis square off
title(sprintf('FZP, on-axis, z0 = %g', 2*z1))
subplot(2,2,3); plot(s, R1(:,128)/254, s, R2(:,128)/254); grid
title('1-D plots of FZP through center'); xlabel('x (a.u.)')
legend(sprintf('z0 = %g', z1), sprintf('z0 = %g', 2*z1), 'Location', 'northeast')
subplot(2,2,4); image(R3); axis square off
title(sprintf('FZP center moves to (x0,y0)=(%g,%g), same depth', x0, -x0))
