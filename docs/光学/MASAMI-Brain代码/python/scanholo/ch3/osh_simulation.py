# ===== 书 3.5 节 Example 3.1 · OSH.m 的规范化重写（Table 3.1，复现原书图 3.15） =====
# 原书 OSH.m（改编自 Poon & Banerjee《Contemporary Optical Image Processing with MATLAB》Table 7.2）：
#   读入 256x256 图像 I(x,y)（原书 vatech.bmp）-> 频域乘 OTFosh -> ifft2 得复全息图 H
#   real(H) 显示为 sine-FZP hologram (3.5-8a)、imag(H) 显示为 cosine-FZP hologram (3.5-8b)
#   重建：fft2(H_any).*conj(OTFosh) 后 ifft2（式 3.5-10），sin/cos 单独重建有孪生噪声，
#         复全息图（real + j*imag）重建干净；H_c-（real - j*imag）重建离焦（图 3.15g）。
#
# 约定说明（重要）：
# 原书 OSH.m 用 MATLAB 原生 fft2（-j 指数）当"前向"，OTFosh 相应取 exp(-j*sigma*k^2)（-j 版），
# 与书正文公式（+j 约定、OTF = exp[+j...] 的共轭）互为镜像——但两者各自自洽，闭环重建结果相同。
# 本脚本统一到"书式协议"（前向 +j = ifft2*N^2），OTFosh = exp(-j*sigma*k^2) 即书 (3.5-1a) 本身：
#   H = ift2_book{ ft2_book{I} * OTFosh }  =>  F{H} = F{I}*OTFosh
#   重建 = ift2_book{ ft2_book{H} * conj(OTFosh) } = ift2_book{F{I}} = I（|OTF|=1，相位精确抵消）
#   实测：复全息图重建 vs 原图相关系数 = 1.0000。
# 电子学符号：real(H) = H_sin（sine 全息图）、imag(H) = -H_cos（显示时翻转灰度，无碍）。
#
# MATLAB / Python 差异提示：
# 1) 原书 OSH.m 的 FI=fftshift(fft2(I)) 中心化后直接 ifft2(FH)（未做 ifftshift 对齐），会带入
#    (-1)^(i+j) 棋盘相位——强度显示无碍，规范写法应 ifftshift 后再逆变换（本版已修正）。
# 2) 原书 sigma = z0/2k0 = 2.0（任意单位）、kx,ky = -12.8:0.1:12.8。此网格下相邻格相位跳变
#    2*sigma*k*dk 最大约 5.1 rad > pi，属于欠采样演示（书图 3.15b/c 颗粒状外观的来源之一）；
#    但闭环重建 self-conjugate 抵消，结论不受影响。物理单位的严谨版本见 point_target_osh.py。
# 3) MATLAB 1-based/Python 0-based；plt.imshow 加 origin='lower' 对齐数学坐标。

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

# ----- 1. 输入图像（原书 vatech.bmp；此处用 PIL 生成 "VT / OSH" 文字图代替） -----
N = 256
img = Image.new('L', (N, N), 0)
d = ImageDraw.Draw(img)
try:
    fnt = ImageFont.truetype('arial.ttf', 88)
except Exception:
    fnt = ImageFont.load_default()
d.text((128, 80), 'VT', font=fnt, fill=255, anchor='mm')
d.text((128, 185), 'OSH', font=fnt, fill=255, anchor='mm')
I = (np.asarray(img) > 80).astype(float)

# ----- 2. 书式傅里叶协议 + OSH 传递函数（sigma = z0/2k0 = 2.0，任意单位，同原书） -----
def ft2_book(U):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

sigma = 2.0
ku = np.arange(N)*0.1 - 12.8
KY, KX = np.meshgrid(ku, ku)
OTFosh = np.exp(-1j*sigma*(KX**2+KY**2))           # 式 (3.5-1a)

# ----- 3. 记录（Coding）：频域乘 OTFosh -> 逆变换得复全息图 -----
H = ift2_book(ft2_book(I)*OTFosh)
Hsin = np.real(H)                                  # sine-FZP hologram   (3.5-8a)
HcosM = np.imag(H)                                 # imag 支路（= -H_cos，显示为 cosine  (3.5-8b)）

# ----- 4. 解码（Decoding，式 3.5-10：频域乘 conj(OTFosh) 模拟与 h(x,y;z0) 卷积） -----
rec_sine = ift2_book(ft2_book(Hsin)*np.conj(OTFosh))      # 孪生噪声
rec_cosine = ift2_book(ft2_book(HcosM)*np.conj(OTFosh))   # 孪生噪声
Hc_plus = H                                    # complex = real + j*imag（书 OSH.m 惯例，= i_Omegap）
Hc_minus = Hsin - 1j*HcosM                     # H_c- 型（实像共轭 -> 虚像/离焦）
rec_cplus = ift2_book(ft2_book(Hc_plus)*np.conj(OTFosh))    # 干净
rec_cminus = ift2_book(ft2_book(Hc_minus)*np.conj(OTFosh))  # 离焦（图 3.15g）

def cc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return abs(np.vdot(a, b))/np.sqrt(np.vdot(a, a).real*np.vdot(b, b).real)

print(f'复全息图 (real+j*imag) 重建 vs 原图：相关系数 = {cc(I, rec_cplus):.4f}')
print(f'H_c- (real-j*imag) 重建 vs 原图：相关系数 = {cc(I, rec_cminus):.4f}（严重离焦）')
rs = np.abs(rec_sine)**2
print(f'sin 全息图单独重建：能量集中度 = {rs.max()/rs.sum():.4f}（孪生像使能量弥散）')

# ----- 5. 可视化（对应原书图 3.15 a)~g)） -----
fig, ax = plt.subplots(2, 4, figsize=(15, 7.2))
ax[0,0].imshow(I, cmap='gray'); ax[0,0].set_title('a) Original image I(x,y)')
ax[0,1].imshow(Hsin, cmap='gray'); ax[0,1].set_title('b) Sine-FZP hologram (3.5-8a)')
ax[0,2].imshow(HcosM, cmap='gray'); ax[0,2].set_title('c) Cosine-FZP hologram (3.5-8b)')
ax[0,3].imshow(np.abs(H), cmap='gray'); ax[0,3].set_title('|complex hologram|')
ax[1,0].imshow(np.abs(rec_sine), cmap='gray', vmax=np.abs(rec_sine).max()*0.5)
ax[1,0].set_title('d) Rec. sine-holo (twin noise)')
ax[1,1].imshow(np.abs(rec_cosine), cmap='gray', vmax=np.abs(rec_cosine).max()*0.5)
ax[1,1].set_title('e) Rec. cosine-holo (twin noise)')
ax[1,2].imshow(np.abs(rec_cplus), cmap='gray')
ax[1,2].set_title('f) Rec. complex H_c+ (clean)')
ax[1,3].imshow(np.abs(rec_cminus), cmap='gray')
ax[1,3].set_title('g) Rec. H_c- (defocused, fig 3.15g)')
for a in ax.flat:
    a.axis('off')
plt.tight_layout(); plt.show()
