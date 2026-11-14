% ===== 书 3.5 节 Example 3.1 · OSH.m 的规范化重写（Table 3.1，复现原书图 3.15） =====
% 原书 OSH.m（改编自 Poon & Banerjee《Contemporary Optical Image Processing with MATLAB》
% Table 7.2, Pages 222-223）：读入 256x256 图像 I(x,y)（原书 vatech.bmp）-> 频域乘 OTFosh
% -> ifft2 得复全息图 H；real(H) 显示为 sine-FZP hologram (3.5-8a)、imag(H) 为 cosine-FZP (3.5-8b)；
% 重建：fft2(H_any).*conj(OTFosh) 后 ifft2（式 3.5-10）。sin/cos 单独重建有孪生噪声；
% 复全息图（real + j*imag）重建干净；H_c-（real - j*imag）重建离焦（图 3.15g）。
%
% 差异提示：
% 1) 原书 OSH.m 用 MATLAB 原生 fft2（-j 指数）当"前向"，OTFosh 相应取 exp(-j*sigma*k^2)，
%    与书正文公式（+j 约定）互为镜像——但各自自洽，闭环重建结果相同。本版统一到书式协议：
%    前向 ifft2*N^2（+j），OTFosh = exp(-j*sigma*k^2) 即书 (3.5-1a) 本身。
%    重建 = ift2{ft2{H}.*conj(OTFosh)} = ift2{F{I}} = I（|OTF|=1，相位精确抵消）。
% 2) 原书 FI=fftshift(fft2(I)) 中心化后直接 ifft2(FH)（未做 ifftshift 对齐），会带入
%    (-1)^(i+j) 棋盘相位——强度显示无碍，规范写法应 ifftshift 后再逆变换（本版已修正）。
% 3) 原书 sigma = z0/2k0 = 2.0（任意单位）、kx,ky = -12.8:0.1:12.8。此网格下相邻格相位跳变
%    2*sigma*k*dk 最大约 5.1 rad > pi，属欠采样演示（书图 3.15b/c 颗粒状外观来源之一）；
%    闭环重建 self-conjugate 抵消，结论不受影响。
% 4) 输入图像：原书读 vatech.bmp；本版若无该文件则用内置几何图形代替（方块+圆环+十字）。

clear
N = 256;
% ---- 输入图像 I(x,y)：优先读 vatech.bmp，否则生成几何测试图 ----
if exist('vatech.bmp', 'file')
    I = double(imread('vatech.bmp'));
    I = I(:,:,1);
    I = I/max(I(:));
else
    [Yg, Xg] = meshgrid((0:N-1)-N/2, (0:N-1)-N/2);
    I = zeros(N);
    I(abs(Xg) < 40 & abs(Yg) < 40) = 1;                    % 中心方块
    I(Xg.^2+Yg.^2 > 45^2 & Xg.^2+Yg.^2 < 60^2) = 1;        % 圆环
    I(abs(Xg) < 3 & abs(Yg) < 100) = 1;                    % 竖线
    I(abs(Yg) < 3 & abs(Xg) < 100) = 1;                    % 横线
end

% ---- 书式傅里叶协议 + OSH 传递函数（sigma = z0/2k0 = 2.0，任意单位，同原书） ----
sigma = 2.0;
ku = (0:N-1)*0.1 - 12.8;
[KY, KX] = meshgrid(ku, ku);
OTFosh = exp(-1j*sigma*(KX.^2+KY.^2));         % 式 (3.5-1a)

% ---- 记录（Coding）：频域乘 OTFosh -> 书式逆变换得复全息图 ----
H = ift2_book(ft2_book(I).*OTFosh);
Hsin = real(H);                                % sine-FZP hologram   (3.5-8a)
HcosM = imag(H);                               % imag 支路（= -H_cos，显示为 cosine  (3.5-8b)）

% ---- 解码（Decoding，式 3.5-10：频域乘 conj(OTFosh) 模拟与 h(x,y;z0) 卷积） ----
rec_sine   = ift2_book(ft2_book(Hsin).*conj(OTFosh));      % 孪生噪声
rec_cosine = ift2_book(ft2_book(HcosM).*conj(OTFosh));     % 孪生噪声
Hc_plus  = H;                                              % complex = real + j*imag（书 OSH.m 惯例）
Hc_minus = Hsin - 1j*HcosM;                                % H_c- 型（实像共轭 -> 虚像/离焦）
rec_cplus  = ift2_book(ft2_book(Hc_plus).*conj(OTFosh));   % 干净
rec_cminus = ift2_book(ft2_book(Hc_minus).*conj(OTFosh));  % 离焦（图 3.15g）

a = I - mean(I(:)); b = rec_cplus - mean(rec_cplus(:));
corr_cplus = abs(sum(a(:).*conj(b(:))))/sqrt(sum(abs(a).^2,'all')*sum(abs(b).^2,'all'));
a = I - mean(I(:)); b = rec_cminus - mean(rec_cminus(:));
corr_cminus = abs(sum(a(:).*conj(b(:))))/sqrt(sum(abs(a).^2,'all')*sum(abs(b).^2,'all'));
rs = abs(rec_sine).^2;
fprintf('复全息图 (real+j*imag) 重建 vs 原图：相关系数 = %.4f\n', corr_cplus);
fprintf('H_c- (real-j*imag) 重建 vs 原图：相关系数 = %.4f（严重离焦）\n', corr_cminus);
fprintf('sin 全息图单独重建：能量集中度 = %.4f（孪生像使能量弥散）\n', max(rs(:))/sum(rs(:)));

% ---- 可视化（对应原书图 3.15 a)~g)） ----
figure(1)
subplot(2,4,1); imshow(I, []); title('a) Original image I(x,y)')
subplot(2,4,2); imshow(real(H), []); title('b) Sine-FZP hologram (3.5-8a)')
subplot(2,4,3); imshow(HcosM, []); title('c) Cosine-FZP hologram (3.5-8b)')
subplot(2,4,4); imshow(abs(H), []); title('|complex hologram|')
subplot(2,4,5); imshow(abs(rec_sine), []); title('d) Rec. sine-holo (twin noise)')
subplot(2,4,6); imshow(abs(rec_cosine), []); title('e) Rec. cosine-holo (twin noise)')
subplot(2,4,7); imshow(abs(rec_cplus), []); title('f) Rec. complex H_c+ (clean)')
subplot(2,4,8); imshow(abs(rec_cminus), []); title('g) Rec. H_c- (defocused)')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function F = ft2_book(U)
% 书式前向变换（+j 指数）
N = size(U, 1);
F = fftshift(ifft2(ifftshift(U)))*N^2;
end

function U = ift2_book(A)
% 书式逆变换（-j 指数）
N = size(A, 1);
U = fftshift(fft2(ifftshift(A)))/N^2;
end
