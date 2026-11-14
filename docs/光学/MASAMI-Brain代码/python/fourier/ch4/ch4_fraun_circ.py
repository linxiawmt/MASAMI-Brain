# ===== 书 4.5 节：fraun_circ —— 圆孔夫琅禾费衍射（Airy 图样） =====
# 书 4.5 节建了两个 M 文件：jinc.m 与 fraun_circ.m。Python 侧不需要单独的
# jinc 文件——公共模块 base.jinc(rho, R=1.0) = J1(2*pi*rho)/rho 与书的 jinc
# 定义完全一致，本章用它包一层 jinc_book 并说明等价性。
# 参数（书 4.5 节）：w = 1 mm（圆孔半径）、lambda = 0.633 um（He-Ne 激光）、
# z = 50 m、观察窗 L = 5 x 1.22*lambda*z/w ≈ 0.193 m（书取 0.2 m）、M = 250。
# 式 (4.35)：I2(x,y) = (w^2/(lambda*z))^2 * jinc(w*rho/(lambda*z))^2。
# MATLAB/Python 差异：
#   - 书的 nthroot(I2,3) 对应 Python 的 I2**(1/3)（3 次方根拉出圆环）；
#   - MATLAB besselj(1,u) 对应 scipy.special.j1(u)（第一类 1 阶贝塞尔）；
#   - imagesc 默认 y 向下 + axis xy 翻回，对应 imshow(origin='lower')。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import coords, jinc   # 公共积木（ch1 起）：jinc(rho, R=1) = 书 ch4 的 jinc


def jinc_book(u):
    """书 4.5 节的 jinc：J1(2*pi*u)/u，u=0 处极限 pi（= base.jinc(u, R=1.0)）。

    书原话："not all jinc functions are the same"——本定义的 jinc(0) = pi，
    是"半径 1 圆盘频谱"约定；base.jinc(rho, R) 的峰值为圆盘面积 pi*R^2，
    R=1 时两者重合（jinc 第一个零点都在 u = 0.61，对应 J1 零点 1.22*pi）。
    """
    return jinc(u, R=1.0)


# ---------- 参数与网格（书 fraun_circ.m） ----------
L, M = 0.2, 250            # 观察窗边长 (m)、样本数
dx = L / M                 # 采样间隔
x, _ = coords(L, M)        # 式 (2.21)：-0.1 : dx : 0.1-dx（coords 的第二个返回值是频率轴！）
y = x                      # y 与 x 同网格（书代码 y=x）
X, Y = np.meshgrid(x, y)   # 'xy' 约定：X 的行 = x 的副本（与 MATLAB 一致）

w = 1e-3          # 圆孔半径 (m)（书注释写 "x half-width"：圆孔的半宽 = 半径）
lam = 0.633e-6    # 波长 (m)，He-Ne
z = 50.0          # 传播距离 (m)
lz = lam * z      # 组合量 lambda*z（书代码变量 lz）

# ---------- 夫琅禾费辐照度（式 4.35，直接用解析 jinc 形式） ----------
I2 = (w ** 2 / lz) ** 2 * jinc_book(w / lz * np.sqrt(X ** 2 + Y ** 2)) ** 2

# ---------- 数值检查 ----------
peak_an = (w ** 2 / lz) ** 2 * np.pi ** 2     # jinc(0) = pi
row = M // 2                                  # y=0 剖面（MATLAB 的 I2(M/2+1,:)）
prof = I2[row, :]
x0_an = 0.61 * lz / w                         # 式 (4.36)：2pi*(w/lz)x = 1.22pi
i0 = int(np.argmin(np.abs(x - x0_an)))
print(f"峰值 I2(0,0) = {I2[row, row]:.6e}；解析 (w^2/(lam*z))^2 * pi^2 = {peak_an:.6e}")
print(f"第一零点解析 0.61*lam*z/w = {x0_an:.5f} m；"
      f"最近网格样本 x = {x[i0]:.5f} m，I2 = {prof[i0]:.2e}（≈0）")
print(f"D_lobe = 1.22*lam*z/w = {1.22 * lz / w:.5f} m；"
      f"L = 5*D_lobe = {5 * 1.22 * lz / w:.4f} m（书据此取 0.2 m）")
print(f"NF = w^2/(lam*z) = {w ** 2 / lz:.4f} << 1"
      f"（书判据 NF<0.1 -> z > 10*w^2/lam = {10 * w ** 2 / lam:.2f} m，取 z=50 m）")
print(f"J1 第一个零点 = 1.22*pi（精确 1.2197*pi）；jinc 第一个零点 u = 0.61")

# ---------- 画图：复现书图 4.3 ----------
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))

axes[0].imshow(I2 ** (1 / 3), extent=[x[0], x[-1], y[0], y[-1]],
               origin="lower", cmap="gray")               # 书：nthroot(I2,3)
axes[0].set_xlabel("x (m)"); axes[0].set_ylabel("y (m)")
axes[0].set_title(r"$I_2^{1/3}$ image: Airy pattern (Fig. 4.3a)")

axes[1].plot(x, prof)                                     # 书：plot(x, I2(M/2+1,:))
axes[1].set_xlabel("x (m)"); axes[1].set_ylabel("Irradiance")
axes[1].set_title("x-axis profile (Fig. 4.3b)")
axes[1].grid(True, alpha=0.3)
# 第一零点参考线（式 4.36/4.37：半宽 0.61*lam*z/w）
for sgn in (-1, 1):
    axes[1].axvline(sgn * x0_an, ls=":", lw=1, color="crimson")
axes[1].annotate(r"$\pm 0.61\lambda z/w$", (x0_an, prof.max() * 0.6),
                 color="crimson", fontsize=9, ha="left")

fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_fraun_circ.png"), dpi=150)
print("已保存 out_fraun_circ.png")
