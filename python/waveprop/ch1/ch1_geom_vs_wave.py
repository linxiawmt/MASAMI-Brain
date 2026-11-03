# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 1 章 =====
几何光学 vs 波动光学：不透明圆盘（水滴模型）挡光后的光强分布。
几何光学给出"绝对阴影"；波动光学（第 1 章的菲涅耳衍射积分，式 1.57）给出
阴影边缘的衍射条纹，以及阴影中心残留的亮斑（泊松斑）。
重写自用户旧稿 几何光学和波动光学对比.py：旧版用 1-(D/λz)·sinc 的启发式公式示意，
本版改为按式 (1.57) 的传递函数法真实传播，物理量纲正确。
用法：直接运行弹窗显示；`python ch1_geom_vs_wave.py <输出目录>` 则保存 PNG。
依赖：numpy + matplotlib。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ===== 1. 参数 =====
N = 1024
dx = 2e-6            # 2 um 像素
wvl = 532e-9
k = 2 * np.pi / wvl
z = 5e-3             # 传播 5 mm（传递函数法窗口内：z < 4 dx^2/λ ≈ 30 mm）

x = (np.arange(N) - N // 2 + 0.5) * dx
X, Y = np.meshgrid(x, x)
R = np.sqrt(X ** 2 + Y ** 2)

D = 100e-6           # 圆盘（水滴）直径 100 um

# ===== 2. 两个模型 =====
# 几何光学：圆盘内光强为 0，其余为 1（绝对阴影）
I_geom = (R > D / 2).astype(float)

# 波动光学：U1 = 1 - circ（平面波减圆盘），用式(1.57)的 TF 形式传播
U1 = (R > D / 2).astype(float)
fx = (np.arange(N) - N // 2) / (N * dx)
FX, FY = np.meshgrid(fx, fx)
H = np.exp(1j * k * z) * np.exp(-1j * np.pi * wvl * z * (FX ** 2 + FY ** 2))
U_ft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(U1)))
U2 = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(U_ft * H)))
I_wave = np.abs(U2) ** 2

# ===== 3. 数值验证 =====
print('===== 几何 vs 波动数值验证 =====')
ctr = N // 2
# (1) 泊松斑：几何光学预言盘心 I=0；波动光学预言轴上振幅恢复到入射值（I=1）。
#     网格离散下泊松斑尺度 sqrt(λz)=1.6 um 小于像素，中心值是它的像素平均残影。
print(f'(1) 盘心光强：几何 = 0，波动数值 = {I_wave[ctr, ctr]:.3f}（理想连续值 1）')
# (2) 阴影边缘：几何光学的跳变 vs 波动光学的条纹（中心行剖面）
row_g, row_w = I_geom[ctr, :], I_wave[ctr, :]
out_sh = (np.abs(x) > D / 2) & (np.abs(x) < D / 2 + 50e-6)   # 阴影外 50 um 环带
print(f'(2) 阴影外 50 um 环带内 波动光强范围 = [{row_w[out_sh].min():.3f}, '
      f'{row_w[out_sh].max():.3f}]（几何恒为 1，衍射产生明暗条纹）')
# (3) 能量守恒：TF 传播前后总能量（|U|^2 求和）应基本不变
print(f'(3) 传播前 sum|U1|^2 = {np.sum(U1**2):.1f}，传播后 sum|U2|^2 = {np.sum(I_wave):.1f}')

# ===== 4. 绘图 =====
def make_figures():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    im = ax1.imshow(I_wave, extent=[x[0] * 1e6, x[-1] * 1e6] * 2, cmap='inferno',
                    vmin=0, vmax=1, origin='lower')
    fig.colorbar(im, ax=ax1, shrink=0.85).set_label('相对光强')
    ax1.add_patch(plt.Circle((0, 0), D / 2 * 1e6, fill=False, color='w', ls='--', lw=0.8))
    ax1.set_xlim(-200, 200); ax1.set_ylim(-200, 200)
    ax1.set_title(f'波动光学：圆盘后 z={z*1e3:.0f} mm 的光强分布')
    ax1.set_xlabel('x (um)'); ax1.set_ylabel('y (um)')

    ax2.plot(x * 1e6, row_g, 'r--', lw=1.5, label='几何光学（绝对阴影）')
    ax2.plot(x * 1e6, row_w, 'b-', lw=1.5, label='波动光学（菲涅耳衍射）')
    ax2.axvspan(-D / 2 * 1e6, D / 2 * 1e6, color='gray', alpha=0.12)
    ax2.set_xlim(-150, 150); ax2.set_ylim(-0.05, 1.15)
    ax2.set_xlabel('x (um)'); ax2.set_ylabel('相对光强')
    ax2.set_title('中心行剖面：两种模型的对比')
    ax2.legend(); ax2.grid(alpha=0.3)

    fig.tight_layout()
    return fig


if __name__ == '__main__':
    if len(sys.argv) > 1:
        outdir = sys.argv[1]
        os.makedirs(outdir, exist_ok=True)
        fig = make_figures()
        fig.savefig(os.path.join(outdir, 'geom_vs_wave.png'), dpi=150)
        print('已保存：geom_vs_wave.png ->', outdir)
    else:
        make_figures()
        plt.show()
