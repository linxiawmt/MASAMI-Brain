% ===== 书 Example 4.2 + Table 4.2 Cryptography.m · 光学扫描加密（OSC）规范化重写 =====
% 加密级（图 4.15 上半）：p2 = δ（针孔，@ω0+Ω），p1 = 加密密钥（@ω0）
%   OTF (4.3-5)：OTF_enc(kx,ky;z^c) = exp[-j z^c(kx²+ky²)/(2k0)] · p1*(-f/k0·kx, -f/k0·ky)
%   加密图像 (4.3-11)：i_enc = F^{-1}{ F{|Γ0|²·exp[j2πr(x,y)]} · exp[-j z^c k²/(2k0)] · exp[-j2πM(-f/k0 k)] }
%   ——双随机相位编码（double-random phase encoding）[Refregier and Javidi (1995)]：
%     输入面随机相位 exp[j2πr] + 频域（光瞳面）随机相位 exp[-j2πM(-u)]，r、M 相互独立。
% 解密级：p1 = δ（针孔，@ω0），p2 = 解密密钥 exp[j2πM(-x,-y)]（@ω0+Ω），扫描针孔物体 @z^d
%   解密信息 (4.3-12)：i_dec = F^{-1}{ exp[+j z^d k²/(2k0)] · exp[+j2πM(-f/k0 k)] }
% DDU（图 4.16）：输出 ∝ F^{-1}{ F{i_enc} × F{i_dec} }（两个复谱直接相乘，无共轭）
%   = |Γ0|²·exp[j2πr] · exp[j(z^d-z^c)k²/(2k0)]·[p1*(-u)p2(u)]
%   条件 1) z^d = z^c（编码距离 = 解码距离）  2) p1*(-x,-y)·p2(x,y) = 1
%   最后取模 |·| 去掉输入面随机相位 exp[j2πr]（|e^{j2πr}|=1）-> 干净文档。
%
% 三种情形（对应书图 4.18；Python 实跑结论）：
%   a) 密钥匹配 + z^d = z^c       -> 相关系数 1.0000（干净）
%   b) 密钥不匹配（猜随机相位）   -> 0.0004（噪声）
%   c) 密钥匹配 + z^d = 1.5 z^c   -> 0.2693（散焦）
%
% 约定说明：原书 Cryptography.m 用 MATLAB 原生 fft2（-j 指数）作"前向"——与本笔记书式协议
% （前向 +j = ifft2*N²，同 ch2/ch3）互为镜像；但 OSC 的解密只依赖 |OTF·P| = 1 的模守恒闭环，
% 两种约定结果相同（延续 ch3 的"解析传递函数敏感/闭环安全"判据）。本版统一书式协议。
%
% 与原书 m-file 的差异：
% 1) 原书读 575×577 的 vt.bmp 后 interp2 到 256×256；本版无该文件时生成几何文字图。
% 2) 原书三重 for 循环生成 OTF；本版向量化（meshgrid 一次生成）。

clear
N = 256;
Dl = 0.02;                        % 物理视场（原书 Dl=0.02）
dx = Dl/N;
zc = 0.30;                        % 编码距离 (m)
lam = 0.6e-6;                     % 波长 (m)
sigma = zc*lam/(4*pi);            % = z^c/(2k0)（同原书 sigma=(zc*ld)/(4*pi)）

x = ((0:N-1) - N/2) * dx;
[X, Y] = meshgrid(x, x);
dfk = 2*pi/(N*dx);
kx1d = ((0:N-1) - N/2) * dfk;
[KY, KX] = meshgrid(kx1d, kx1d);

% ---- 待加密文档（原书 vt.bmp "Virginia Tech"；无该文件则生成几何测试图） ----
if exist('vt.bmp', 'file')
    C = imread('vt.bmp');
    C = double(C(:,:,1));
    [xs, ys] = meshgrid(linspace(1, size(C,2), N), linspace(1, size(C,1), N));
    I0 = interp2(1:size(C,2), 1:size(C,1), C, xs, ys);
    I0 = I0 / max(I0(:));
else
    I0 = zeros(N);
    I0(abs(X) < 0.3e-3 & abs(Y) < 0.7e-3 & abs(Y) > 0.35e-3) = 1;   % 上下两条横杠
    I0(abs(X + 0.55e-3) < 0.12e-3 & abs(Y) < 0.75e-3) = 1;          % V 左斜
    I0(abs(X + 0.2e-3)  < 0.12e-3 & abs(Y) < 0.75e-3) = 1;          % V 右斜
    I0(abs(X - 0.35e-3) < 0.03e-3 & abs(Y) < 0.75e-3) = 1;          % T 竖
    I0(abs(X - 0.35e-3) < 0.55e-3 & abs(Y - 0.6e-3) < 0.05e-3) = 1; % T 横
end

% ---- 双随机相位编码（式 4.3-11） ----
rng(7);                                    % 固定随机种子便于复现（原书 rand 不固定）
R = rand(N);                               % 输入面随机相位 r(x,y) ~ U(0,1)
Mkey = rand(N);                            % 加密密钥 M(x,y)（p1 = exp[j2πM]，频域取 M(-u)）
RPM = exp(1j*2*pi*R);                      % exp[j2πr]
P = exp(-1j*2*pi*Mkey);                    % 频域密钥 = exp[-j2πM(-u)]（随机矩阵，等价）
OTF_enc = exp(-1j*sigma*(KX.^2+KY.^2));    % exp[-j z^c k²/(2k0)]（式 4.3-5 的相位项）

R1 = I0 .* RPM;                            % 文档 × 输入面随机相位
E = ft2_book(R1) .* OTF_enc .* P;          % 加密频谱（式 4.3-11 花括号内）
e = ift2_book(E);                          % 加密图像 i(x,y;z^c)（密文，发送）

% ---- 解密级信息（式 4.3-12）与 DDU ----
OTF_dec   = exp(1j*sigma*(KX.^2+KY.^2));            % exp[+j z^d k²/(2k0)]，z^d = z^c
OTF_dec15 = exp(1j*1.5*sigma*(KX.^2+KY.^2));        % z^d = 1.5 z^c
Key_info = conj(OTF_enc .* P);                      % F{i_dec}（p2 = P*，满足条件 2）
M3 = rand(N);
Key_mis = conj(OTF_enc .* exp(-1j*2*pi*M3));        % 猜错密钥（独立随机相位）

de     = ift2_book(E .* Key_info);                  % a) 匹配密钥 + z^d = z^c
de_mis = ift2_book(E .* Key_mis);                   % b) 密钥不匹配
de_15  = ift2_book(E .* conj(OTF_dec15 .* P));      % c) 匹配密钥 + z^d = 1.5 z^c

fprintf('--- 解密质量（|DDU输出| vs 原文档 I0 的相关系数） ---\n');
fprintf('a) 匹配密钥 + z^d = z^c      : %.4f（干净）\n', cc(abs(de), I0));
fprintf('b) 密钥不匹配 + z^d = z^c    : %.4f（噪声）\n', cc(abs(de_mis), I0));
fprintf('c) 匹配密钥 + z^d = 1.5·z^c  : %.4f（散焦）\n', cc(abs(de_15), I0));
fprintf('条件核对：|OTF_enc·P| = %.6f（模恒 1 -> DDU 乘积中密钥项精确对消）\n', mean(abs(OTF_enc.*P), 'all'));
fprintf('密文 |e| 与文档的相关系数 = %.4f（密文强度不含文档信息 -> 窃听无效）\n', cc(abs(e), I0));

% ---- 可视化（对应书图 4.17/4.18） ----
figure(1)
subplot(2,4,1); imagesc(x*1e3, x*1e3, I0); axis image off; colormap gray
title('a) Original document |G_0|^2')
subplot(2,4,2); imagesc(x*1e3, x*1e3, real(R1)); axis image off
title('b) Re{doc × exp(j2πr)}')
subplot(2,4,3); imagesc(x*1e3, x*1e3, imag(R1)); axis image off
title('c) Im{doc × exp(j2πr)}')
subplot(2,4,4); imagesc(x*1e3, x*1e3, abs(e)); axis image off
title('d) |encrypted| (sent via wireless)')
subplot(2,4,5); imagesc(x*1e3, x*1e3, abs(de)); axis image off
title(sprintf('a) matched key, zd=zc (cc=%.3f)', cc(abs(de), I0)))
subplot(2,4,6); imagesc(x*1e3, x*1e3, abs(de_mis)); axis image off
title(sprintf('b) mismatched key (cc=%.3f)', cc(abs(de_mis), I0)))
subplot(2,4,7); imagesc(x*1e3, x*1e3, abs(de_15)); axis image off
title(sprintf('c) zd=1.5zc (cc=%.3f)', cc(abs(de_15), I0)))
subplot(2,4,8); imagesc(x*1e3, x*1e3, angle(E)); axis image off
title('spectrum phase of encrypted image')

% ---- 局部函数（脚本内定义，须置于文件末尾，R2016b+ 支持） ----
function F = ft2_book(U)
% 书式前向变换（+j 指数）
n = size(U, 1);
F = fftshift(ifft2(ifftshift(U))) * n^2;
end

function u = ift2_book(A)
% 书式逆变换（-j 指数）
n = size(A, 1);
u = fftshift(fft2(ifftshift(A))) / n^2;
end

function r = cc(a, b)
% 去均值复相关系数的模
a = a - mean(a(:));  b = b - mean(b(:));
r = abs(sum(a(:).*conj(b(:)))) / sqrt(sum(abs(a).^2, 'all') * sum(abs(b).^2, 'all'));
end
