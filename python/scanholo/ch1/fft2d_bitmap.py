# ===== 书 1.1 例 1.2 · Table 1.4 fft2Dbitmap_image.m 的 Python 版 =====
# 对一幅 256x256 位图做二维傅里叶变换并显示幅度谱。
# 原书读 triangle.bmp（画图板随手画的三角形）；为了让你零依赖跑通，
# 这里用 numpy 现场画一个同样的三角形代替 imread。
#
# MATLAB / Python 差异提示：
# 1) np.fft.fft2 与 MATLAB fft2 同为"负指数"约定，可直接照抄；
#    fftshift 也同名同义（零频搬到中心）。
# 2) 原书 image(10*(abs(256*FI))) 是把幅度放大后按 0~255 灰度显示；
#    matplotlib 中等价写法是 imshow 归一化后的幅度谱。

import numpy as np
import matplotlib.pyplot as plt

# ----- 生成 256x256 三角形位图（替代 I=imread('triangle.bmp','bmp')） -----
N = 256
I = np.zeros((N, N))
yy, xx = np.mgrid[0:N, 0:N]
inside = (yy > 160) & (yy < 160 + (xx - 128) * 0.9) & (xx > 48) & (xx < 208)
I[inside] = 1.0                      # 白色三角形，黑色背景（幅度 1 = 白）

# ----- 书 figure(1)：显示输入图像 -----
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].imshow(I, cmap='gray', origin='lower')
ax[0].set_title('a) Bitmap image of a triangle')
ax[0].axis('off')

# ----- 计算傅里叶变换并中心化 -----
FI = np.fft.fft2(I)                  # I=I(:,:,1)：取一个通道（这里本来就是单通道）
FI = np.fft.fftshift(FI)
scale = 1.0 / np.abs(FI).max()       # max2=max(max1); scale=1.0/max2
FI = FI * scale

# ----- 书 figure(2)：显示变换的绝对值（对数拉伸便于观看星芒） -----
ax[1].imshow(20 * np.log10(np.abs(FI) + 1e-6), cmap='gray', origin='lower')
ax[1].set_title('b) |Transform| of a) (dB display)')
ax[1].axis('off')

plt.tight_layout()
plt.show()
