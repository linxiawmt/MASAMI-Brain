% ===== 书 3.5 节 · 点目标扫描全息：sin/cos FZP 全息图与复全息图重建（式 3.5-1 ~ 3.5-10） =====
% 物理设置：p1(x,y)=1（光瞳全开 -> L1 后焦面上点源，频率 w0）、p2(x,y)=delta（-> 平面波，w0+Omega）。
% 点目标 |Gamma_0|^2 = delta(x,y)delta(z-z0) 位于距扫描镜 z0 处。
% 主结果：
%   H_sin(x,y) = (k0/2*pi*z0) sin[(k0/2*z0)(x^2+y^2)]           (3.5-4a)
%   H_cos(x,y) = (k0/2*pi*z0) cos[(k0/2*z0)(x^2+y^2)]           (3.5-4b)
%   H_c+ = H_cos + j*H_sin = j*i_Omegap  ->  与 h(x,y;z0) 卷积聚焦成点（conc = 1.0000）
%   H_c- = H_cos - j*H_sin = -j*conj(i_Omegap) -> 同样卷积得到 2*z0 处菲涅耳衍射图（corr = 1.0000）
%   sin 全息图单独重建：主像只占一半能量（另一半给孪生像/离焦斑，实测 0.501）——孪生像问题。
%
% 差异提示：
% 1) 全书最大的坑：书的前向变换是 +j 指数，而 MATLAB 的 fft2 是 -j 指数。
%    严格按书公式实现前向变换必须写 ifft2(U)*N^2；逆变换写 fft2(A)/N^2（局部函数在文件末尾）。
% 2) e^{-j sigma k^2} 不是带限函数：离散记录的 sin/cos 全息图含 Poisson 旁瓣混叠
%    （书 fig 3.15b/c 颗粒状外观的数值根源），但闭环重建用 conj(OTFosh) 精确抵消，不受影响。
% 3) MATLAB 1-based：中心像素为 (257,257)。

clear
lam = 0.6328e-6;                 % He-Ne 红光波长（m）
k0 = 2*pi/lam;                   % 波数
z0 = 0.30;                       % 点目标到扫描镜的距离（m）
N = 512;                         % 网格
dx = 10e-6;                      % 像素尺寸（m）
sigma = z0/(2*k0);               % = z0/(2*k0)，式 (3.5-1a) 的二次相位系数

k = 2*pi*((0:N-1)-N/2)/(N*dx);   % 中心格角频率
[KY, KX] = meshgrid(k, k);

% ---- OSH 传递函数 (3.5-1a) 与点目标记录 ----
OTFosh = exp(-1j*sigma*(KX.^2+KY.^2));
I_pt = zeros(N); I_pt(N/2+1, N/2+1) = 1;       % 点目标（中心 1 像素）
iOmegap = ift2_book(ft2_book(I_pt).*OTFosh);   % i_Omegap = F^{-1}{F{|Gamma0|^2} OTF_osh}  (3.4-15)

Hsin = real(iOmegap);                          % H_sin = Re{i_Omegap}   (3.5-3a)
Hcos = -imag(iOmegap);                         % H_cos = -Im{i_Omegap}  (3.5-3b，电子学负号)
Hc_plus  = Hcos + 1j*Hsin;                     % H_c+ = H_cos + j*H_sin = j*i_Omegap    (3.5-5)
Hc_minus = Hcos - 1j*Hsin;                     % H_c- = H_cos - j*H_sin = -j*conj(i_Omegap)

% ---- 数字重建（式 3.5-10：与 h(x,y;z0) 卷积 = 频域乘 OTF_osh*） ----
rec_cp = ift2_book(ft2_book(Hc_plus).*conj(OTFosh));     % 聚焦 -> 点（实像）
rec_cm = ift2_book(ft2_book(Hc_minus).*conj(OTFosh));    % 离焦 -> 2*z0 菲涅耳衍射
rec_sp = ift2_book(ft2_book(Hsin).*conj(OTFosh));        % sin 全息图单独 -> 点 + 孪生

rcp = abs(rec_cp).^2;
[~, pk] = max(rcp(:));  [pk_r, pk_c] = ind2sub(size(rcp), pk);
conc = rcp(pk)/sum(rcp(:));
rsm = abs(rec_sp).^2;
ratio_sin = rsm(pk)/sum(rsm(:));

% ---- 对照：H_c- 重建 vs 原物在 2*z0 处的菲涅耳衍射（书 p.88 文字结论） ----
H2 = exp(-1j*k0*2*z0).*exp(1j*(KX.^2+KY.^2)*(2*z0)/(2*k0));
fres2 = ift2_book(ft2_book(I_pt).*H2);
a = abs(rec_cm); b = abs(fres2);
a = a - mean(a(:)); b = b - mean(b(:));
corr_cm = abs(sum(a(:).*conj(b(:))))/sqrt(sum(abs(a).^2,'all')*sum(abs(b).^2,'all'));

fprintf('H_c+ 重建：峰值位置 (%d,%d)（点目标在 (%d,%d)），能量集中度 = %.4f\n', ...
        pk_r, pk_c, N/2+1, N/2+1, conc);
fprintf('sin 全息图单独重建：主像能量占比 = %.3f（另一半进入孪生像）\n', ratio_sin);
fprintf('|H_c- 重建| vs |2*z0 菲涅耳衍射|：相关系数 = %.4f\n', corr_cm);

% ---- 可视化 ----
Xg = ((0:N-1)-N/2)*dx;
[Yg, Xg] = meshgrid(Xg, Xg);
x = Xg(1,:)*1e3;                                 % mm
figure(1)
subplot(2,3,1); imagesc(x, x, sin(k0*(Xg.^2+Yg.^2)/(2*z0))); axis image off; colormap gray
title('sin-FZP of a point (Eq. 3.5-4a, analytic)')
subplot(2,3,2); imagesc(x, x, cos(k0*(Xg.^2+Yg.^2)/(2*z0))); axis image off
title('cos-FZP of a point (Eq. 3.5-4b, analytic)')
subplot(2,3,3); imagesc(x, x, rcp.^0.3); axis image off
title(sprintf('H_c+ * h(z0): focused (conc = %.2f)', conc))
subplot(2,3,4); imagesc(x, x, rsm.^0.3); axis image off
title('sin-FZP alone: point + twin (50/50)')
subplot(2,3,5); imagesc(x, x, abs(rec_cm).^0.5); axis image off
title('H_c- * h(z0): defocused')
subplot(2,3,6); imagesc(x, x, abs(fres2).^0.5); axis image off
title('Fresnel diffraction @2z0 (match)')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function F = ft2_book(U)
% 书式前向变换（+j 指数）：F(kx,ky) = 双重积分 U exp[j(kx x+ky y)] dx dy
N = size(U, 1);
F = fftshift(ifft2(ifftshift(U)))*N^2;
end

function U = ift2_book(A)
% 书式逆变换（-j 指数）：U(x,y) = (1/4pi^2) 双重积分 F exp[-j(kx x+ky y)] dkx dky
N = size(A, 1);
U = fftshift(fft2(ifftshift(A)))/N^2;
end
