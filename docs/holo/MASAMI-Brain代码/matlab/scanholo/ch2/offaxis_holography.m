% ===== 书 2.5.2 · 离轴点源全息：记录（式 2.5-7）与再现（式 2.5-8） =====
% 改编自你已有的 holography.m（规范化重写 + 统一到书约定）：
% 点源位于 (0,0,z0)（胶片在 z=0，点源在 z=-z0），离轴平面参考波以记录角 theta 入射。
% 记录：t(x,y) = |psi_r + psi_o|^2 = A + B sin[ (k0/2z0)(x^2+y^2) + k0 x sin(theta) ]
% 再现：用同一参考波照射全息图，菲涅耳传播：
%     零级波沿 theta 方向继续出射                       (2.5-8a)
%     实像（孪生像）在胶片前方 z=+z0，横向 x=-2 z0 sin(theta)   (2.5-8b)
%     虚像（原始点源位置）在胶片后方 z=-z0，横向 (0,0)          (2.5-8c)
%
% 差异提示：
% 1) 本脚本全程严格书约定：前向变换 ifft2*N^2、逆变换 fft2/N^2（函数在文件末尾）。
%    你旧版 holography.m 用 fft2 配 exp(+j k^2 z/2k0)：fft2 与书约定互为镜像，两块
%    强度图依然正确（相当于结果左右镜像），但复振幅是镜像的，且实/虚像标注要与镜像
%    约定配套才对——本版用书协议后，实像严格落在 x=-2 z0 sin(theta)（式 2.5-8b）。
% 2) 物波按书 (2.5-7) 取 psi_o = (j k0/2 pi z0) exp[-j k0 r^2/(2 z0)]（书同样省略
%    全局相位 exp(-j k0 z0)）；参考波振幅 a 取 |psi_o| 最大值，使条纹对比度最大。
% 3) 采样检查：载频 fc = sin(theta)/lambda0，本例约 158 cycle/cm，dx≈9.8 um，
%    每条条纹约 6.5 个采样点，满足 Nyquist（呼应 2.5.3 的分辨率讨论）。

clear
L = 1;                  % 胶片/探测器尺寸：1 cm x 1 cm
N = 1024;               % 采样点数（1024 便于解析载频条纹）
dx = L/(N-1);
X = (0:N-1)*dx - L/2;
Y = (0:N-1)*dx - L/2;

lambda0 = 0.6328e-4;    % He-Ne 红光波长（cm）
k0 = 2*pi/lambda0;      % 波数
z0 = 10;                % 点源到胶片的距离（cm）
theta = 0.01;           % 记录角（rad）；2*z0*sin(theta)=0.2 cm 在窗口内

% ---- 频域网格（DFT 真实间距） ----
dk = 2*pi/(N*dx);
kx = (-N/2:N/2-1)*dk;
ky = (-N/2:N/2-1)*dk;
[KXX, KYY] = meshgrid(kx, ky);

% ---- 记录（Coding）：物波球面波 + 离轴参考平面波 ----
[XX, YY] = meshgrid(X, Y);
r2 = XX.^2 + YY.^2;
psi_o = (1j*k0/(2*pi*z0))*exp(-1j*k0*r2/(2*z0));    % 物波（发散球面波，式 2.5-7）
a = max(abs(psi_o(:)));                             % 参考波振幅（1:1 条纹对比度）
psi_r = a*exp(1j*k0*XX*sin(theta));                 % 离轴参考波（书式 +j）

Hologram = abs(psi_r + psi_o).^2;                   % 全息图 t(x,y)
A0 = a^2 + (k0/(2*pi*z0))^2;                        % 式 (2.5-7) 的 A
B0 = 2*a*(k0/(2*pi*z0));                            % 式 (2.5-7) 的 B = a k0/(pi z0)

% ---- 再现（Decoding）：参考波照射 + 菲涅耳传播到 +z0 / -z0 ----
U_trans   = Hologram .* psi_r;                      % psi_rc * t(x,y)
U_real    = prop_book(U_trans, +z0, k0, KXX, KYY, N);   % 传播 +z0：实像（孪生像）平面
U_virtual = prop_book(U_trans, -z0, k0, KXX, KYY, N);   % 反向传播 -z0：虚像平面
I_real    = abs(U_real).^2;
I_virtual = abs(U_virtual).^2;

x_real_theory = -2*z0*sin(theta);                   % 式 (2.5-8b)：delta(x+2 z0 sinθ, y)
[~, ir] = max(I_real(:));    [rr_, cr_] = ind2sub(size(I_real), ir);
[~, iv] = max(I_virtual(:)); [rv_, cv_] = ind2sub(size(I_virtual), iv);
fprintf('实像峰值位置  x = %+.4f cm （理论 %+.4f cm, y=0）\n', X(cr_), x_real_theory);
fprintf('虚像峰值位置  x = %+.4f cm （理论 %+.4f cm, y=0）\n', X(cv_), 0);

% ---- 可视化 ----
cx = abs(X) < 0.06;                                 % 全息图中心 1.2mm 放大看载频条纹
figure(1)
subplot(2,2,1); imagesc(X(cx), Y(cx), Hologram(cx, cx)); colormap gray; axis image
title('Off-axis hologram t(x,y), zoomed center'); xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,2,2); imagesc(X, Y, I_real.^0.35); colormap gray; axis image
title(sprintf('Real image (twin) at z=+%g cm, x=%+.3f cm', z0, X(cr_)));
xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,2,3); imagesc(X, Y, I_virtual.^0.35); colormap gray; axis image
title(sprintf('Virtual image at z=-%g cm, x=%+.3f cm', z0, X(cv_)));
xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,2,4); plot(X, I_real(round(N/2),:)/max(I_real(:)), 'r-'); hold on
plot(X, I_virtual(round(N/2),:)/max(I_virtual(:)), 'b--');
xline(x_real_theory, 'k:');
xlim([-0.4 0.4]); grid; legend('real-image plane (z=+z0)', 'virtual-image plane (z=-z0)')
title('1-D intensity slices (y=0)'); xlabel('x (cm)')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function U = prop_book(U0, z, k0, KXX, KYY, N)
% 角谱传播（式 2.3-5 精确形式，自动滤消逝波）：z 可正可负
root = 1 - (KXX.^2 + KYY.^2)/k0^2;
root(root < 0) = 0;
H = exp(-1j*k0*z*sqrt(root));
A = fftshift(ifft2(ifftshift(U0)))*N^2;   % 书式前向变换（+j）
U = fftshift(fft2(ifftshift(A.*H)))/N^2;  % 书式逆变换（-j）
end
