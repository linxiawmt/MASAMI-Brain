# ===== 书 1.2 例 1.6 · Table 1.5 correlation.m 的 Python 版 =====
# 频域实现二维相关（图案识别）：
#     f(x,y) (+) h(x,y) = IFFT{ F*(kx,ky) H(kx,ky) }        [书式 (1.2-7)(1.2-8)]
# 原书读两幅 smiley.bmp 做自相关（图 1.14 中心亮斑 = 两图"匹配"）。
# 这里用 numpy 现场画笑脸替代 imread；把 I2 换成别的图案即可复现图 1.15
# 的互相关（中心无亮斑）。
#
# MATLAB / Python 差异提示：
# 1) np.fft.fft2 / ifft2 / fftshift 与 MATLAB 同名函数约定一致，可逐行照抄。
# 2) 共轭放在哪个谱上（FI1.*conj(FI2) 还是 conj(FI1).*FI2）只决定输出是否
#    镜像翻转；对实图像取 abs 显示后看不出差别。数值验证见 verify_ch1.py。

import numpy as np
import matplotlib.pyplot as plt

# ----- 现场生成 256x256 笑脸位图（替代 I=imread('smiley.bmp','bmp')） -----
def make_smiley():
    N = 256
    yy, xx = np.mgrid[0:N, 0:N]                # 行 yy、列 xx（0..255）
    r = np.sqrt((xx - 128.0)**2 + (yy - 128.0)**2)
    img = np.zeros((N, N))
    img[np.abs(r - 90) < 4] = 1.0                                  # 外圈脸廓
    for ex in (93, 163):                                           # 两只眼睛
        img[(xx - ex)**2 + (yy - 88)**2 < 9**2] = 1.0
    ang = np.arctan2(yy - 118.0, xx - 128.0)                       # 嘴巴（下方圆弧）
    rm = np.sqrt((xx - 128.0)**2 + (yy - 118.0)**2)
    img[(np.abs(rm - 55) < 4) & (ang < -0.35) & (ang > -2.8)] = 1.0
    return img

# ----- 输入图像 1（参考图像）与输入图像 2（待识别图像） -----
I1 = make_smiley()        # I1=imread('smiley.bmp','bmp'); I1=I1(:,:,1)
I2 = make_smiley()        # 原书两幅同名文件 => 自相关；换成别的图案 => 互相关

fig, ax = plt.subplots(1, 3, figsize=(13, 4.6))
ax[0].imshow(I1, cmap='gray', origin='lower'); ax[0].set_title('Input image 1 (reference)'); ax[0].axis('off')
ax[1].imshow(I2, cmap='gray', origin='lower'); ax[1].set_title('Input image 2 (to recognize)'); ax[1].axis('off')

# ----- 频域计算相关 -----
FI1 = np.fft.fft2(I1)
FI1 = FI1 / np.abs(FI1).max()          # max2=max(max1); scale=1.0/max2; FI1=FI1.*scale
FI2 = np.fft.fft2(I2)
FI2 = FI2 / np.abs(FI2).max()

FPR = FI1 * np.conj(FI2)               # calculating correlation（书中 FPR=FI1.*conj(FI2)）
PR = np.fft.ifft2(FPR)
PR = np.fft.fftshift(PR)
PR = np.abs(PR) / np.abs(PR).max()     # 归一化后取模显示

ax[2].imshow(PR, cmap='gray', origin='lower')
ax[2].set_title('Correlation result (bright spot = match)')
ax[2].axis('off')

pk = np.unravel_index(np.argmax(PR), PR.shape)
print(f'相关峰位置 (row, col) = {pk}，应落在图像中心 (128, 128)')
plt.tight_layout()
plt.show()
