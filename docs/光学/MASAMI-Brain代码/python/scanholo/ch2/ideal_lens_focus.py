# ===== 书 2.4.1 · 理想透镜相位变换 + 角谱传播：动态聚焦仿真 =====
# 改编自你已有的 Ideal Lens.m / Ideal Lens.py（规范化重写）：
# 输入图案（同心圆靶盘+十字线）紧贴理想透镜（d0=0），用式 (2.3-13) 的
# 菲涅耳传递函数传播到 z = 0.5f, f, 1.5f，观察聚焦过程；
# 并对比 d0=0（透镜贴合物，焦面场带相位弯曲，式 2.4-3）与
# d0=f（物放前焦面，焦面场相位平坦的严格傅里叶变换，式 2.4-6）。
#
# 核心公式：
#     透镜相位变换  t_f(x,y) = exp[ +j k0 (x^2+y^2) / (2f) ]        (2.4-1)
#     菲涅耳传递函数 H(kx,ky;z) = exp(-j k0 z) exp[ j (kx^2+ky^2) z / (2 k0) ]  (2.3-13)
#
# 采样警告（旧版踩过的坑）：透镜啁啾相位的局部空间频率为 k0*x/f，随 x 线性增长；
# 它必须小于网格 Nyquist 频率 pi/dx，即 |x| < f*lambda/(2 dx)。旧参数 L=1cm、
# N=256（dx≈3.9um）只允许 |x|<0.08cm，窗外全部混叠，焦面场被污染。
# 本版取 L=0.4cm、N=1024（dx≈3.9e-4 cm），允许 |x|<0.81cm，图案半径 0.17cm，
# 有 4 倍余量。这与 2.5.3 的"记录介质分辨率"是同一件事：相位越陡，采样要越密。
#
# MATLAB / Python 差异提示：
# 1) 书的前向变换是 +j 指数，严格实现用 ifft2*N^2（fft2 是 -j 指数，见 asm_propagation.py
#    注释 1）。你旧版 Ideal Lens.m 用 fft2 + 书式传播子 exp(+j k^2 z/2k0)：fft2 与书约定
#    互为镜像，强度图仍然正确（相当于结果左右镜像），但复振幅/相位是镜像的——本版
#    统一改成书协议，相位穹顶对比才是严格对的。
# 2) meshgrid 默认 indexing='xy'，广播 r2 = X[:,None]**2 + Y[None,:]**2 等价双层 for。
# 3) 频率格间距用 DFT 真实值 dk = 2*pi/(N*dx)。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 物理参数 -----
L = 0.4                    # 空间尺寸：0.4 cm x 0.4 cm（保证透镜啁啾被充分采样，见上）
N = 1024                   # 网格大小
dx = L / (N - 1)
X = (np.arange(N)) * dx - L / 2
Y = (np.arange(N)) * dx - L / 2

lam = 0.6328e-4            # 红光波长：0.6328 um（cm）
k0 = 2 * np.pi / lam       # 波数
f = 10.0                   # 透镜焦距：10 cm

# ----- 2. 频域网格标定（角谱算法，DFT 真实间距） -----
dk = 2 * np.pi / (N * dx)
kx = (np.arange(N) - N // 2) * dk
ky = (np.arange(N) - N // 2) * dk
KX, KY = np.meshgrid(kx, ky)

# ----- 3. 制造输入图案（同心圆靶盘与十字线，单位 cm） -----
R = np.zeros((N, N))
XX, YY = np.meshgrid(X, Y)
rr = np.sqrt(XX**2 + YY**2)
R[(rr > 0.030) & (rr < 0.050)] = 1
R[(rr > 0.085) & (rr < 0.105)] = 1
R[(rr > 0.140) & (rr < 0.160)] = 1
R[np.abs(XX) < 1.5 * dx] = 1               # 十字线
R[np.abs(YY) < 1.5 * dx] = 1
R[rr > 0.170] = 0                          # 裁剪边缘

# ----- 4. 书式傅里叶协议封装（前向 +j = ifft2*N^2，逆向 -j = fft2/N^2） -----
def ft2_book(U):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

def prop_book(U0, z):
    """精确角谱传播（式 2.3-5）：先 FT，频域乘 H，再 IFT"""
    root = 1 - (KX**2 + KY**2) / k0**2
    root[root < 0] = 0
    H = np.exp(-1j * k0 * z * np.sqrt(root))
    return ift2_book(ft2_book(U0) * H)

# ----- 5. 紧贴透镜后的光场 = 物体乘透镜相位 t_f（式 2.4-1） -----
r2 = X[:, None]**2 + Y[None, :]**2
U_behind_lens = R * np.exp(1j * (k0 / (2 * f)) * r2)

# ----- 6. 传播到 z = 0.5f, f, 1.5f 观察聚焦 -----
U_z = [prop_book(U_behind_lens, zc) for zc in (0.5 * f, f, 1.5 * f)]

# ----- 7. d0=f 对照组：物放前焦面 {[t*h(x,y;d0=f)]*t_f}*h(x,y;f)（式 2.4-4，取 d0=f） -----
U_df = prop_book(prop_book(R, f) * np.exp(1j * (k0 / (2 * f)) * r2), f)

# ----- 8. 可视化 -----
fig, ax = plt.subplots(2, 3, figsize=(13, 8))
xms = [X[0], X[-1]]

ax[0, 0].imshow(R, extent=[*xms, *xms], cmap='gray', origin='lower')
ax[0, 0].set_title('Input target pattern R(x,y)')

panels = [(0, 1, 'Intensity z = 0.5f (log)', U_z[0]),
          (0, 2, 'Intensity z = f (focus, log)', U_z[1]),
          (1, 0, 'Intensity z = 1.5f (log)', U_z[2])]
for rr_, cc_, title, U in panels:
    I = np.abs(U)**2; I /= I.max()
    ax[rr_, cc_].imshow(np.log10(I * 100 + 1), extent=[*xms, *xms], cmap='hot', origin='lower')
    ax[rr_, cc_].set_title(title)

# 相位对比：d0=0 焦面场带相位弯曲（式 2.4-3）vs d0=f 相位平坦（式 2.4-6）
mask0 = np.abs(U_z[1])**2 > 1e-4 * np.abs(U_z[1]).max()**2
maskf = np.abs(U_df)**2 > 1e-4 * np.abs(U_df).max()**2
ph0 = np.where(mask0, np.angle(U_z[1]), np.nan)
phf = np.where(maskf, np.angle(U_df), np.nan)
ax[1, 1].imshow(ph0, extent=[*xms, *xms], cmap='jet', origin='lower')
ax[1, 1].set_title('Phase at focus, d0 = 0 (curved, Eq.2.4-3)')
ax[1, 2].imshow(phf, extent=[*xms, *xms], cmap='jet', origin='lower')
ax[1, 2].set_title('Phase at focus, d0 = f (flat, Eq.2.4-6)')

plt.tight_layout()
plt.show()
