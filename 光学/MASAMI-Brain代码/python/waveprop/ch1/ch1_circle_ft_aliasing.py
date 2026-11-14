# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 1 章 =====
圆孔的二维傅里叶变换 = 夫琅禾费衍射（艾里斑），以及粗采样带来的频谱混叠。
重写自用户旧稿 Fourier_Transfer_2D.m / circle_digital.m / 二维傅里叶变换.py / 圆孔离散化.py：
  - cv2 画圆改为解析 circ(R<=D/2)，去掉 cv2 依赖；
  - 细/粗两套网格合并到一个脚本，补上与连续理论的定量核对。
用法：直接运行弹窗显示；`python ch1_circle_ft_aliasing.py <输出目录>` 则保存 PNG。
依赖：numpy + matplotlib。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (注册 3d projection)

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ===== 1. 参数 =====
D = 100e-6          # 圆孔直径 100 um（对应原书 circ 函数，Appendix A）
wvl = 532e-9        # 波长 532 nm

# --- A. 细网格：逼近连续圆孔 ---
N_f = 512
dx_f = 1e-6                       # 1 um 像素，圆孔直径占 100 像素
x_f = (np.arange(N_f) - N_f // 2 + 0.5) * dx_f
X_f, Y_f = np.meshgrid(x_f, x_f)
U_f = (np.sqrt(X_f ** 2 + Y_f ** 2) <= D / 2).astype(float)   # circ(2r/D)

# --- B. 粗网格：同一个圆孔，像素放大 4 倍 ---
N_c = 128
dx_c = 4e-6                       # 4 um 像素，圆孔直径只占 25 像素 -> 边缘锯齿
x_c = (np.arange(N_c) - N_c // 2 + 0.5) * dx_c
X_c, Y_c = np.meshgrid(x_c, x_c)
U_c = (np.sqrt(X_c ** 2 + Y_c ** 2) <= D / 2).astype(float)

# ===== 2. 居中二维傅里叶变换（ft2，含面积元 dx^2）=====
def ft2(g, delta):
    """G(fx,fy) = FFT2(g) 的中心化排列 * delta^2，近似连续傅里叶变换。
    坐标约定 x = (n - N/2 + 1/2) delta（胞元中心），频率轴 fx = (k - N/2)/(N delta)。
    MATLAB 同名实现见笔记，注意 fftshift 在 MATLAB 里对每维同样适用。"""
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g))) * delta ** 2

G_f = ft2(U_f, dx_f)
I_f = np.abs(G_f) ** 2
I_f /= I_f.max()
f_f = (np.arange(N_f) - N_f // 2) / (N_f * dx_f)          # 频率轴 (1/m)

G_c = ft2(U_c, dx_c)
I_c = np.abs(G_c) ** 2
I_c /= I_c.max()
f_c = (np.arange(N_c) - N_c // 2) / (N_c * dx_c)

# ===== 3. 数值验证 =====
print('===== 圆孔傅里叶变换数值验证 =====')
# (1) 直流值：FT{circ(2r/D)}(0) = 圆孔面积 = pi D^2 / 4
dc_th = np.pi * D ** 2 / 4
dc_num = np.abs(G_f[N_f // 2, N_f // 2])
print(f'(1) 直流分量理论 pi*D^2/4 = {dc_th:.4e}，数值 = {dc_num:.4e}，相对偏差 {abs(dc_num-dc_th)/dc_th:.2e}')

# (2) 艾里斑第一暗环：J1 零点 3.8317 -> f = 3.8317/(pi*D)；三点抛物线细化极小位置
f1_th = 3.8317 / (np.pi * D)
row = I_f[N_f // 2, :]                       # 沿 fx 轴（fy=0）的剖面
half = N_f // 2 + 1                          # f>0 半轴
rings = [i for i in range(half, N_f - 1) if row[i] < row[i - 1] and row[i] < row[i + 1]]
i0 = rings[0]
df_bin = f_f[1] - f_f[0]
den = row[i0 - 1] - 2 * row[i0] + row[i0 + 1]
f1_num = f_f[i0] + 0.5 * df_bin * (row[i0 - 1] - row[i0 + 1]) / (den + 1e-30)
print(f'(2) 第一暗环理论 f = {f1_th/1e3:.3f} mm^-1，数值（细化）= {f1_num/1e3:.3f} mm^-1，'
      f'相对偏差 {abs(f1_num-f1_th)/f1_th:.2e}')

# (3) 粗网格的两类失真（沿 fx 剖面对比细网格）：
#     像素化边缘把光散射进艾里暗环（低频段暗环被"填亮"）；
#     频谱以 1/dx_c 为周期复制，在近奈奎斯特端（|f|>100 mm^-1）出现混叠尾巴。
fx_c_mm = f_c * 1e-3
I_f_interp = np.interp(fx_c_mm, f_f * 1e-3, row / row[half])   # 细网格归一化剖面
prof_c = I_c[N_c // 2, :] / I_c[N_c // 2, N_c // 2]
i_ring = int(np.argmin(np.abs(fx_c_mm - f1_th * 1e-3)))
fill = prof_c[i_ring] / max(I_f_interp[i_ring], 1e-12)
high = np.abs(fx_c_mm) > 100
err_high = np.max(np.abs(prof_c[high] - I_f_interp[high]))
print(f'(3) 第一暗环处 粗/细强度比 = {fill:.1f}（像素化把暗环填亮，细网格处接近 0）；'
      f'|f|>100 mm^-1 频谱复制尾巴最大绝对偏差 = {err_high:.2e}')

# ===== 4. 绘图 =====
def make_figures():
    fig = plt.figure(figsize=(15, 9.5))

    ax1 = fig.add_subplot(2, 3, 1)
    ax1.imshow(U_f, extent=[x_f[0] * 1e6, x_f[-1] * 1e6] * 2, cmap='gray', origin='lower')
    ax1.set_title('(a) 细网格圆孔（连续近似）')
    ax1.set_xlabel('x (um)'); ax1.set_ylabel('y (um)')

    ax2 = fig.add_subplot(2, 3, 2)
    im2 = ax2.imshow(np.log10(I_f + 1e-6), extent=[f_f[0] * 1e-3, f_f[-1] * 1e-3] * 2,
                     cmap='gray', vmin=-5, vmax=0, origin='lower')
    fig.colorbar(im2, ax=ax2, shrink=0.85).set_label('log$_{10}$(I/I$_{max}$)')
    ax2.set_title('(b) 艾里斑频谱（对数尺度）')
    ax2.set_xlabel('$f_x$ (mm$^{-1}$)'); ax2.set_ylabel('$f_y$ (mm$^{-1}$)')

    ax3 = fig.add_subplot(2, 3, 3, projection='3d')
    s = 4                                     # 3D 面抽稀，仅中心 1/4 区域
    c = slice(N_f // 2 - 32 * s, N_f // 2 + 32 * s, s)
    FXg, FYg = np.meshgrid(f_f[c] * 1e-3, f_f[c] * 1e-3)
    ax3.plot_surface(FXg, FYg, I_f[c, c], cmap='jet', edgecolor='none')
    ax3.set_title('(c) 艾里斑 3D 光强')
    ax3.set_xlabel('$f_x$ (mm$^{-1}$)'); ax3.set_ylabel('$f_y$ (mm$^{-1}$)')

    ax4 = fig.add_subplot(2, 3, 4)
    ax4.imshow(U_c, extent=[x_c[0] * 1e6, x_c[-1] * 1e6] * 2, cmap='gray', origin='lower')
    ax4.set_title('(d) 粗网格圆孔（边缘锯齿）')
    ax4.set_xlabel('x (um)'); ax4.set_ylabel('y (um)')

    ax5 = fig.add_subplot(2, 3, 5)
    im5 = ax5.imshow(np.log10(I_c + 1e-6), extent=[f_c[0] * 1e-3, f_c[-1] * 1e-3] * 2,
                     cmap='gray', vmin=-5, vmax=0, origin='lower')
    fig.colorbar(im5, ax=ax5, shrink=0.85).set_label('log$_{10}$(I/I$_{max}$)')
    ax5.set_title('(e) 粗网格频谱：高频端出现混叠')
    ax5.set_xlabel('$f_x$ (mm$^{-1}$)'); ax5.set_ylabel('$f_y$ (mm$^{-1}$)')

    ax6 = fig.add_subplot(2, 3, 6)
    ax6.plot(f_f * 1e-3, row / row[half], 'b-', lw=1.2, label='细网格 (dx=1 um)')
    ax6.plot(fx_c_mm, prof_c, 'r.-', lw=0.8, ms=4, label='粗网格 (dx=4 um)')
    ax6.axvline(f1_th * 1e-3, color='g', ls=':', lw=1.2,
                label='理论第一暗环 12.2 mm$^{-1}$')
    ax6.set_xlim(-60, 60); ax6.set_ylim(1e-5, 1.1)
    ax6.set_yscale('log')
    ax6.set_xlabel('$f_x$ (mm$^{-1}$)'); ax6.set_ylabel('归一化光强')
    ax6.set_title('沿 $f_x$ 剖面对比（对数纵轴）')
    ax6.legend(); ax6.grid(alpha=0.3, which='both')

    fig.tight_layout()
    return fig


if __name__ == '__main__':
    if len(sys.argv) > 1:
        outdir = sys.argv[1]
        os.makedirs(outdir, exist_ok=True)
        fig = make_figures()
        fig.savefig(os.path.join(outdir, 'circle_ft_aliasing.png'), dpi=150)
        print('已保存：circle_ft_aliasing.png ->', outdir)
    else:
        make_figures()
        plt.show()
