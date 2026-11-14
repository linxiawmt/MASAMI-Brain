# ===== 书 3.7 节：二维矩形函数、meshgrid、图像显示与二维 FFT =====
# 参数：wx = 0.1 m、wy = 0.05 m 的矩形（0.2 × 0.1 m），L = 2 m、M = 200/边
# -> dx = 0.01 m，200×200 阵列对应 2 × 2 m。
# MATLAB/Python 差异：
#   - np.meshgrid 默认 indexing='xy' 与 MATLAB meshgrid 一致（X 的行是 x
#     的副本、Y 的列是 y 的副本）；写 'ij' 两维会互换；
#   - MATLAB g(M/2+1,:) 选第 M/2+1 行 = y=0 的行（行 <-> y！）；Python
#     0-based 是 g[M//2, :]；
#   - plt.imshow 默认 y 向下，origin='lower' 翻转（MATLAB 用 axis xy）；
#   - 二维必须 fft2（np.fft.fft2）：一维 fft 作用于二维阵列只变换一行行
#     向量，得到"条纹"状错误结果（书 3.8 节原话）。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  （注册 3D 投影，新版本可省）

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import coords, sinc   # 公共积木（ch2 起）


def rect_book(x):
    """书 3.1 节的 rect 编码：abs(x) <= 1/2 -> 1（见 ch3_fft_example.py 注释）。"""
    return (np.abs(x) <= 0.5).astype(float)


wx, wy = 0.1, 0.05        # 矩形 x/y 半宽 (m)
L, M = 2.0, 200           # 边长 (m)、每边样本数
dx = L / M                # 采样间隔 (m)
x, _ = coords(L, M)       # x 坐标
y = x                     # y 坐标
X, Y = np.meshgrid(x, y)  # X 的行 = x 的副本；Y 的列 = y 的副本
g = rect_book(X / (2 * wx)) * rect_book(Y / (2 * wy))   # 二维矩形（可分离乘积）

# ---------- 二维 FFT ----------
# 书：g0=fftshift(g); G0=fft2(g0)*dx^2; G=fftshift(G0)——偶数 M 下
# fftshift 与 ifftshift 等价；本文件按第 2 章的严谨配对写 ifftshift。
G0 = np.fft.fft2(np.fft.ifftshift(g)) * dx ** 2   # 二维 FFT × dx^2 定标
G = np.fft.fftshift(G0)                           # 居中
_, fx = coords(L, M)                              # 频率坐标（式 2.22；coords 首个返回值是空域 x！）
fy = fx                                           # 两维相同

# ---------- 数值检查 ----------
row = M // 2                       # y=0 的行（MATLAB g(M/2+1,:)）
# 书 rect 编码（边界取 1）下，x 向 |X|<=0.1 命中 21 个样本、y 向 11 个：
# wx=0.1、wy=0.05 都是"偶数样本宽"，按 3.1 节的约定各多给一个边界样本。
nx = int(rect_book(x / (2 * wx)).sum())   # 21
ny = int(rect_book(y / (2 * wy)).sum())   # 11
# 离散精确谱（Dirichlet 核）：FFT 剖面应与它机器精度重合
with np.errstate(invalid="ignore", divide="ignore"):
    D = np.sin(nx * np.pi * fx * dx) / np.sin(np.pi * fx * dx)
D[fx == 0] = nx                                   # 0/0 极限
slice_exact = dx ** 2 * ny * np.abs(D)
slice_an = 4 * wx * wy * np.abs(sinc(2 * wx * fx))  # 连续矩形解析谱
in_main = np.abs(fx) <= 25
err_exact = np.max(np.abs(np.abs(G[row, :][in_main]) - slice_exact[in_main]))
err_an = np.max(np.abs(np.abs(G[row, :][in_main]) - slice_an[in_main]))
print(f"x/y 向值为 1 的样本数 = {nx} x {ny}（边界取 1 的奇数样本约定）")
print(f"DC 样本 |G(0,0)| = {abs(G[row, row]):.6f} = 离散面积 "
      f"{nx}·{ny}·dx² = {nx * ny * dx ** 2:.4f}")
print(f"  （连续矩形面积 4wx·wy = {4 * wx * wy:.3f}；书图 3.11(b) 峰值 ~0.023 "
      "就是离散面积，")
print("    不是混叠——wx、wy 恰好都是偶数样本宽，rect 各多给一个边界样本）")
print(f"fy=0 剖面 vs 离散精确谱（Dirichlet）最大差异 = {err_exact:.2e}（机器精度）")
print(f"fy=0 剖面 vs 连续解析谱最大差异 = {err_an:.4f}"
      "（DC 处 0.0231 vs 0.020 + 零点错位）")

# ---------- 显示：图像 + 切片 + 曲面 + 频域剖面（书图 3.10/3.11） ----------
fig = plt.figure(figsize=(12.5, 8.4))

ax1 = fig.add_subplot(2, 2, 1)
ax1.imshow(g, extent=[x[0], x[-1], y[0], y[-1]],
           origin="lower", cmap="gray")          # origin='lower' = MATLAB axis xy
ax1.set_xlabel("x (m)"); ax1.set_ylabel("y (m)")
ax1.set_title("g(x,y): 2D rect (Fig. 3.10a)")

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(x, g[row, :], "o-", ms=3, color="crimson", mfc="none")  # y=0 行切片
ax2.set_xlim(-1, 1); ax2.set_ylim(0, 1.2)
ax2.set_xlabel("x (m)")
ax2.set_title("x-slice through g[M//2,:] (Fig. 3.10b)")

ax3 = fig.add_subplot(2, 2, 3, projection="3d")
FXf, FYf = np.meshgrid(fx, fy)                     # 曲面用频率网格
ax3.plot_surface(FXf, FYf, np.abs(G), cmap="gray",
                 rstride=2, cstride=2, linewidth=0)
ax3.set_xlabel("fx (cyc/m)"); ax3.set_ylabel("fy (cyc/m)")
ax3.set_title(r"$|G(f_X,f_Y)|$ surface (Fig. 3.11a)")

ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(fx, np.abs(G[row, :]), label="FFT")
ax4.plot(fx, slice_exact, ":", lw=1, label="discrete exact (Dirichlet)")
ax4.plot(fx, slice_an, "--", lw=1, label=r"$4w_xw_y\,\mathrm{sinc}(2w_x f_x)$")
ax4.set_xlim(-50, 50); ax4.set_xlabel("fx (cyc/m)")
ax4.set_title("fx slice: FFT vs two references (Fig. 3.11b)")
ax4.legend(fontsize=8)

for ax in (ax1, ax2, ax4):
    ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_fft2_example.png"), dpi=150)
print("已保存 out_fft2_example.png")
