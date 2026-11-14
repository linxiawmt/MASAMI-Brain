% ===== 书 2.4.1 · 理想透镜相位变换 + 角谱传播：动态聚焦仿真 =====
% 改编自你已有的 Ideal Lens.m（规范化重写）：
% 输入图案（同心圆靶盘+十字线）紧贴理想透镜（d0=0），用式 (2.3-13) 的
% 菲涅耳传递函数传播到 z = 0.5f, f, 1.5f，观察聚焦过程；
% 并对比 d0=0（透镜贴合物，焦面场带相位弯曲，式 2.4-3）与
% d0=f（物放前焦面，焦面场相位平坦的严格傅里叶变换，式 2.4-6）。
%
% 核心公式：
%     透镜相位变换  t_f(x,y) = exp[ +j k0 (x^2+y^2) / (2f) ]        (2.4-1)
%     菲涅耳传递函数 H(kx,ky;z) = exp(-j k0 z) exp[ j (kx^2+ky^2) z / (2 k0) ]  (2.3-13)
%
% 采样警告（旧版踩过的坑）：透镜啁啾相位的局部空间频率为 k0*x/f，随 x 线性增长；
% 它必须小于网格 Nyquist 频率 pi/dx，即要求 |x| < f*lambda/(2 dx)。旧参数 L=1cm、
% N=256（dx≈3.9um）只允许 |x|<0.08cm，窗外全部混叠，焦面场被污染。
% 本版取 L=0.4cm、N=1024（dx≈3.9e-4 cm），允许 |x|<0.81cm，图案半径 0.17cm，
% 有 4 倍余量。这与 2.5.3 的"记录介质分辨率"是同一件事：相位越陡，采样要越密。
%
% 差异提示：书的前向变换是 +j 指数，严格实现用 ifft2*N^2（fft2 是 -j 指数）。
% 你旧版 Ideal Lens.m 用 fft2 + 书式传播子 exp(+j k^2 z/2k0)：fft2 与书约定互为
% 镜像，强度图仍然正确（相当于结果左右镜像），但复振幅/相位是镜像的——本版统一
% 改成书协议，相位穹顶对比才是严格对的。局部函数在文件末尾（R2016b+）。

clear
L = 0.4;                % 空间尺寸：0.4 cm x 0.4 cm（保证透镜啁啾被充分采样，见上）
N = 1024;               % 网格大小
dx = L/(N-1);
X = (0:N-1)*dx - L/2;
Y = (0:N-1)*dx - L/2;
[XX, YY] = meshgrid(X, Y);

lambda0 = 0.6328e-4;    % 红光波长（cm）
k0 = 2*pi/lambda0;      % 波数
f = 10;                 % 透镜焦距：10 cm

% ---- 频域网格标定（角谱算法，DFT 真实间距） ----
dk = 2*pi/(N*dx);
kx = (-N/2:N/2-1)*dk;
ky = (-N/2:N/2-1)*dk;
[KXX, KYY] = meshgrid(kx, ky);

% ---- 制造输入图案（同心圆靶盘与十字线，单位 cm） ----
rr = sqrt(XX.^2 + YY.^2);
R = zeros(N);
R(rr > 0.030 & rr < 0.050) = 1;
R(rr > 0.085 & rr < 0.105) = 1;
R(rr > 0.140 & rr < 0.160) = 1;
R(abs(XX) < 1.5*dx) = 1;        % 十字线
R(abs(YY) < 1.5*dx) = 1;
R(rr > 0.170) = 0;              % 裁剪边缘

% ---- 书式傅里叶协议 + 传播（函数在文件末尾） ----
r2 = XX.^2 + YY.^2;
tf = exp(1j*(k0/(2*f))*r2);             % 透镜相位变换 t_f（式 2.4-1）

U_behind_lens = R .* tf;                % 紧贴透镜后的光场（d0=0）
z_list = [0.5*f, f, 1.5*f];
U_z = zeros(N, N, 3);
for i = 1:3
    U_z(:,:,i) = prop_book(U_behind_lens, z_list(i), k0, KXX, KYY, N);
end

% d0=f 对照组：物放前焦面 {[t*h(x,y;d0=f)]*t_f}*h(x,y;f)（式 2.4-4，取 d0=f）
U_df = prop_book(prop_book(R, f, k0, KXX, KYY, N).*tf, f, k0, KXX, KYY, N);

% ---- 可视化 ----
figure(1)
subplot(2,3,1); imagesc(X, Y, R); colormap gray; axis image
title('Input target pattern R(x,y)'); xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,3,2); imagesc(X, Y, log10(abs(U_z(:,:,1)).^2/max(abs(U_z(:,:,1)).^2(:))*100+1));
colormap(hot(256)); axis image
title('Intensity z = 0.5f (log)'); xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,3,3); imagesc(X, Y, log10(abs(U_z(:,:,2)).^2/max(abs(U_z(:,:,2)).^2(:))*100+1));
axis image
title('Intensity z = f (focus, log)'); xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,3,4); imagesc(X, Y, log10(abs(U_z(:,:,3)).^2/max(abs(U_z(:,:,3)).^2(:))*100+1));
axis image
title('Intensity z = 1.5f (log)'); xlabel('x (cm)'); ylabel('y (cm)')

% 相位对比：d0=0 焦面场带相位弯曲（式 2.4-3）vs d0=f 相位平坦（式 2.4-6）
I0 = abs(U_z(:,:,2)).^2; If = abs(U_df).^2;
mask0 = I0 > 1e-4*max(I0(:));
maskf = If > 1e-4*max(If(:));
ph0 = angle(U_z(:,:,2)); ph0(~mask0) = NaN;
phf = angle(U_df);       phf(~maskf) = NaN;
subplot(2,3,5); imagesc(X, Y, ph0); colormap jet; axis image
title('Phase at focus, d0 = 0 (curved, Eq.2.4-3)'); xlabel('x (cm)'); ylabel('y (cm)')
subplot(2,3,6); imagesc(X, Y, phf); colormap jet; axis image
title('Phase at focus, d0 = f (flat, Eq.2.4-6)'); xlabel('x (cm)'); ylabel('y (cm)')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function A = ft2_book(U, N)
% 书式前向变换（+j 指数）
A = fftshift(ifft2(ifftshift(U)))*N^2;
end

function U = ift2_book(A, N)
% 书式逆变换（-j 指数）
U = fftshift(fft2(ifftshift(A)))/N^2;
end

function U = prop_book(U0, z, k0, KXX, KYY, N)
% 精确角谱传播（式 2.3-5）：先 FT，频域乘 H，再 IFT
root = 1 - (KXX.^2 + KYY.^2)/k0^2;
root(root < 0) = 0;
H = exp(-1j*k0*z*sqrt(root));
U = ift2_book(ft2_book(U0, N).*H, N);
end
