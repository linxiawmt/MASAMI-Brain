# ===== 拓展实战（编者，原书无）· 同轴数字全息的背景相除法 =====
# 改编自你已有的 holography.m / 孪生像.py / 背景相除法.py（规范化重写，去 OpenCV 依赖）：
# 过冷水滴（振幅 0.95 + 相位 1.2 rad 的圆盘）经 z=5cm 角谱传播生成同轴全息图，
# 相机面叠加高斯激光背景与两粒灰尘 -> 背景相除 I_raw/I_bg -> sqrt 估计复振幅 -> 反向传播重建。
#
# 与本章 OSH 的关系（笔记 3.5 节"孪生像问题"对照）：
# - 背景相除只能去掉乘性背景（灰尘/照明不均），对同轴全息的孪生像虚影无能为力——
#   重建振幅图里水滴周围仍有发散的"孪生像波纹"。
# - OSH 路线（本章）：一次 2-D 扫描同时记录 sin/cos 两张全息图，复数相加 H_cos + j*H_sin
#   即可在同轴记录下彻底消除孪生像（式 3.5-5）。
#
# 严格书约定：前向变换 ifft2*N^2（+j 指数）、逆变换 fft2/N^2（-j 指数）——见 ch2 asm_propagation.py。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 版用 regionprops 做连通域定位；此处用简单阈值 + 质心公式（numpy 足够）。
# 2) imagesc/imshow 坐标原点差异同前几章。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 参数与物体/背景 -----
N = 512
lam = 532e-9
dx = 5e-6
L = N*dx
z_prop = 0.05
k0 = 2*np.pi/lam
x = (np.arange(N)-N//2)*dx
y = (np.arange(N)-N//2)*dx
X, Y = np.meshgrid(x, y)

w0 = 1.5e-3                                        # 高斯激光背景束腰
I_laser = np.exp(-2*(X**2+Y**2)/w0**2)
t_noise = np.ones((N, N))
dust1 = ((X-0.5e-3)**2 + (Y-0.5e-3)**2) <= (20e-6)**2
dust2 = ((X+0.6e-3)**2 + (Y-0.4e-3)**2) <= (15e-6)**2
t_noise[dust1 | dust2] = 0.2                       # 两粒灰尘（透过率 0.2）
I_bg = I_laser*t_noise

U_obj = np.ones((N, N), dtype=complex)
mask = (X**2 + Y**2) <= (50e-6)**2                 # 半径 50 um 的过冷水滴
U_obj[mask] = 0.95*np.exp(1j*1.2)                  # 振幅 0.95、相位 1.2 rad

# ----- 2. 书式傅里叶协议 + 角谱正/反向传播 -----
def ft2_book(U):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

df = 1/L
fx = (np.arange(N)-N//2)*df
FY, FX = np.meshgrid(fx, fx)
KX, KY = 2*np.pi*FX, 2*np.pi*FY
root = 1 - (KX**2+KY**2)/k0**2
root[root < 0] = 0                                 # 滤除消逝波
H_fwd = np.exp(-1j*k0*z_prop*np.sqrt(root))
H_bwd = np.exp(+1j*k0*z_prop*np.sqrt(root))        # 反向传播：z 取 -z_prop

U_det = ift2_book(ft2_book(U_obj)*H_fwd)           # 传播到相机面
I_raw = np.abs(U_det)**2 * I_bg                    # 相机记录（含背景）
I_proc = I_raw/I_bg                                # 背景相除（乘性噪声消除）
U_rec = ift2_book(ft2_book(np.sqrt(I_proc))*H_bwd) # 反向传播重建（sqrt 强度当复振幅估计）

# ----- 3. 简单定位（阈值 + 质心，代替 OpenCV 连通域） -----
amp = np.abs(U_rec)
amp_n = (amp - amp.min())/(amp.max()-amp.min())
binary = amp_n < 0.45
ys, xs = np.nonzero(binary)
if len(xs) > 0:
    cx, cy = xs.mean(), ys.mean()
    d_equiv = 2*dx*np.sqrt(binary.sum()/np.pi)*1e6
    print(f'水滴质心：x = {(cx-N/2)*dx*1e3:+.3f} mm, y = {(cy-N/2)*dx*1e3:+.3f} mm, '
          f'等效直径 = {d_equiv:.1f} um（真值 100 um）')
    print('注意：重建振幅中水滴周围仍可见孪生像虚影波纹——背景相除无法消除孪生像。')

# ----- 4. 可视化 -----
fig, ax = plt.subplots(2, 3, figsize=(14, 8))
ext = [x[0]*1e3, x[-1]*1e3, y[0]*1e3, y[-1]*1e3]
ax[0,0].imshow(np.abs(U_obj), extent=ext, cmap='gray')
ax[0,0].set_title('1. Object amplitude |U_obj| (z=0)')
ax[0,1].imshow(I_raw, extent=ext, cmap='gray')
ax[0,1].set_title('2. Raw hologram (laser bg + dust)')
ax[0,2].imshow(I_proc, extent=ext, cmap='gray')
ax[0,2].set_title('3. Background-divided hologram')
ax[1,0].imshow(amp, extent=ext, cmap='gray')
ax[1,0].set_title('4. Reconstructed amplitude (twin-image haze)')
ax[1,1].imshow(np.angle(U_rec), extent=ext, cmap='jet')
ax[1,1].set_title('5. Reconstructed phase (rad)')
ax[1,2].imshow(amp, extent=ext, cmap='gray')
if len(xs) > 0:
    ax[1,2].plot((cx-N/2)*dx*1e3, (cy-N/2)*dx*1e3, 'ro', ms=10, mfc='none', mew=2)
    ax[1,2].set_title('6. Centroid detection')
for a in ax.flat:
    a.set_xlabel('x (mm)'); a.set_ylabel('y (mm)')
plt.tight_layout(); plt.show()
