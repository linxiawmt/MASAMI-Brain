# ===== 书 Example 4.2 + Table 4.2 Cryptography.m · 光学扫描加密（OSC）规范化重写 =====
# 加密级（图 4.15 上半）：p2 = δ（针孔，@ω0+Ω），p1 = 加密密钥（@ω0）
#   OTF (4.3-5)：OTF_enc(kx,ky;z^c) = exp[-j z^c(kx²+ky²)/(2k0)] · p1*(-f/k0·kx, -f/k0·ky)
#   加密图像 (4.3-11)：i_enc = F^{-1}{ F{|Γ0|²·exp[j2πr(x,y)]} · exp[-j z^c k²/(2k0)] · exp[-j2πM(-f/k0 k)] }
#   ——双随机相位编码（double-random phase encoding）[Refregier and Javidi (1995)]：
#     输入面随机相位 exp[j2πr] + 频域（光瞳面）随机相位 exp[-j2πM(-u)]，r、M 相互独立。
# 解密级：p1 = δ（针孔，@ω0），p2 = 解密密钥 exp[j2πM(-x,-y)]（@ω0+Ω），扫描针孔物体 @z^d
#   解密信息 (4.3-12)：i_dec = F^{-1}{ exp[+j z^d k²/(2k0)] · exp[+j2πM(-f/k0 k)] }
# DDU（图 4.16）：输出 ∝ F^{-1}{ F{i_enc} × F{i_dec} }（两个复谱直接相乘，无共轭）
#   = |Γ0|²·exp[j2πr] · exp[j(z^d-z^c)k²/(2k0)]·[p1*(-u)p2(u)]
#   条件 1) z^d = z^c（编码距离 = 解码距离）  2) p1*(-x,-y)·p2(x,y) = 1
#   最后取模 |·| 去掉输入面随机相位 exp[j2πr]（|e^{j2πr}|=1）-> 干净文档。
#
# 三种情形（对应书图 4.18）：
#   a) 密钥匹配 + z^d = z^c       -> 干净重建（本跑 |de| 与原文档相关系数 ~1.00）
#   b) 密钥不匹配（猜随机相位）   -> 噪声
#   c) 密钥匹配 + z^d = 1.5 z^c   -> 散焦（z^c 与 z^d 是两个独立的密钥自由度）
#
# 约定说明：原书 Cryptography.m 用 MATLAB 原生 fft2（-j 指数）作"前向"——与本笔记书式协议
# （前向 +j = ifft2*N²，同 ch2/ch3）互为镜像；但 OSC 的解密只依赖 |OTF·P| = 1 的模守恒闭环，
# 两种约定结果相同（延续 ch3 的"解析传递函数敏感/闭环安全"判据）。本版统一书式协议。
#
# MATLAB / Python 差异提示：
# 1) 原书 m-file 读 575×577 的 vt.bmp 后 interp2 到 256×256；本版用 PIL 生成 "VT/OSC" 文字图。
# 2) MATLAB 1-based / Python 0-based；imshow 原点差异同前几章。

import sys
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

# ----- 1. 参数（同原书 Cryptography.m：L=1, Dl=0.02, N=255->本版 N=256, zc=0.3, λ=0.6 µm） -----
N = 256
Dl = 0.02                       # 物理视场（原书 Dl=0.02）
dx = Dl / N
zc = 0.30                       # 编码距离 (m)
lam = 0.6e-6                    # 波长 (m)
sigma = zc * lam / (4 * np.pi)  # = z^c/(2k0)（同原书 sigma=(zc*ld)/(4*pi)）

x = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(x, x)
df = 1 / (N * dx)
kx1d = (np.arange(N) - N // 2) * 2 * np.pi * df
KY, KX = np.meshgrid(kx1d, kx1d)

def ft2_book(U):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

# ----- 2. 待加密文档（原书 vt.bmp "Virginia Tech"；本版生成文字图） -----
img = Image.new('L', (N, N), 0)
d = ImageDraw.Draw(img)
try:
    fnt = ImageFont.truetype('arial.ttf', 52)
except Exception:
    fnt = ImageFont.load_default()
d.text((128, 82), 'VT', font=fnt, fill=255, anchor='mm')
d.text((128, 178), 'OSC', font=fnt, fill=255, anchor='mm')
I0 = (np.asarray(img) > 80).astype(float)          # |Γ0(x,y;z^c)|²（0/1 文档）

# ----- 3. 双随机相位编码（式 4.3-11） -----
rng = np.random.default_rng(7)
R = rng.random((N, N))                              # 输入面随机相位 r(x,y) ~ U(0,1)
Mkey = rng.random((N, N))                           # 加密密钥 M(x,y)（p1 = exp[j2πM]，频域取 M(-u)）
RPM = np.exp(1j * 2 * np.pi * R)                    # exp[j2πr]
P = np.exp(-1j * 2 * np.pi * Mkey)                  # 频域密钥 = exp[-j2πM(-u)]（随机矩阵，等价）
OTF_enc = np.exp(-1j * sigma * (KX**2 + KY**2))     # exp[-j z^c k²/(2k0)]   （式 4.3-5 的相位项）

R1 = I0 * RPM                                       # 文档 × 输入面随机相位
E = ft2_book(R1) * OTF_enc * P                      # 加密频谱（式 4.3-11 花括号内）
e = ift2_book(E)                                    # 加密图像 i(x,y;z^c)（密文，发送）

# ----- 4. 解密级信息（式 4.3-12）与 DDU -----
OTF_dec = np.exp(1j * sigma * (KX**2 + KY**2))      # exp[+j z^d k²/(2k0)]，z^d = z^c
OTF_dec15 = np.exp(1j * 1.5 * sigma * (KX**2 + KY**2))   # z^d = 1.5 z^c
Key_info = np.conj(OTF_enc * P)                     # F{i_dec} = OTF_dec·p2((f/k0)k)（p2 = P*，满足条件 2）
M3 = rng.random((N, N))
Key_mis = np.conj(OTF_enc * np.exp(-1j * 2 * np.pi * M3))   # 猜错密钥（独立随机相位）

de = ift2_book(E * Key_info)                        # a) 匹配密钥 + z^d = z^c
de_mis = ift2_book(E * Key_mis)                     # b) 密钥不匹配
de_15 = ift2_book(E * np.conj(OTF_dec15 * P))       # c) 匹配密钥 + z^d = 1.5 z^c

def cc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return abs(np.vdot(a, b)) / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real)

print('--- 解密质量（|DDU输出| vs 原文档 I0 的相关系数） ---')
print(f'a) 匹配密钥 + z^d = z^c      : {cc(np.abs(de), I0):.4f}   （干净）')
print(f'b) 密钥不匹配 + z^d = z^c    : {cc(np.abs(de_mis), I0):.4f}   （噪声）')
print(f'c) 匹配密钥 + z^d = 1.5·z^c  : {cc(np.abs(de_15), I0):.4f}   （散焦）')
print('条件核对：|OTF_enc·P| = '
      f'{np.abs(OTF_enc * P).mean():.6f}（模恒 1 -> DDU 乘积中密钥项精确对消）')
print('密文 |e| 与文档的相关系数 = '
      f'{cc(np.abs(e), I0):.4f}（密文强度不含文档信息 -> 窃听无效）')

# ----- 5. 可视化（对应书图 4.17/4.18） -----
fig, ax = plt.subplots(2, 4, figsize=(15, 7.2))
ax[0, 0].imshow(I0, cmap='gray'); ax[0, 0].set_title('a) Original document |G0|^2')
ax[0, 1].imshow(np.real(R1), cmap='gray'); ax[0, 1].set_title('b) Re{doc x exp(j2πr)}')
ax[0, 2].imshow(np.imag(R1), cmap='gray'); ax[0, 2].set_title('c) Im{doc x exp(j2πr)}')
ax[0, 3].imshow(np.abs(e), cmap='gray'); ax[0, 3].set_title('d) |encrypted| (sent via wireless)')
ax[1, 0].imshow(np.abs(de), cmap='gray'); ax[1, 0].set_title(f"a) matched key, zd=zc (cc={cc(np.abs(de), I0):.3f})")
ax[1, 1].imshow(np.abs(de_mis), cmap='gray'); ax[1, 1].set_title(f"b) mismatched key (cc={cc(np.abs(de_mis), I0):.3f})")
ax[1, 2].imshow(np.abs(de_15), cmap='gray'); ax[1, 2].set_title(f"c) zd=1.5zc (cc={cc(np.abs(de_15), I0):.3f})")
ax[1, 3].imshow(np.angle(E), cmap='twilight', vmin=-np.pi, vmax=np.pi)
ax[1, 3].set_title('spectrum phase of encrypted image')
for a in ax.flat:
    a.axis('off')
plt.tight_layout()
plt.savefig('sim-osc-cryptography.png', dpi=150, bbox_inches='tight')
plt.show()
