# ===== 书 3.9 节习题 3.1（三角函数）与 3.3（圆函数）参考解答 =====
# 3.1：底半宽 0.1 m 的三角函数，L = 2、M = 200。书提示的编码：
#       T=1-abs(x); ask=abs(x)<=1; out=T.*mask;
#     其中 ask 与 mask 是同一个变量（原书笔误，两个名字混用，本文件统一
#     写成 mask）。按 3.1(b) 的尺度（底半宽 w = 0.1）缩放后，
#     解析频谱：tri(x/w) <-> w sinc^2(w fx)，w = 0.1 -> 0.1 sinc^2(0.1 fx)。
# 3.3：半径 0.015 m 的圆盘，L = 0.2、M = 200（dx = 0.001，半径恰 15 个样本）。
#     频谱 = jinc 形：R*J1(2*pi*R*rho)/rho，峰值 = 面积 pi*R^2 ≈ 7.07e-4。
#     （半径约定见第 1 章勘误 1：本解按"半径 = R"的 Goodman 约定，
#      与表 1.3 的 circ 变换对一致。）
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import coords, sinc, tri, circle, jinc   # 公共积木（ch2 起）

# ---------- 习题 3.1：三角函数 ----------
L, M = 2.0, 200
dx = L / M
x, fx = coords(L, M)
w = 0.1                                # 底半宽 (m)
T = 1 - np.abs(x) / w                  # 书提示第 1 行（按习题尺度缩放）
mask = np.abs(x) <= w                  # 书提示第 2 行（原书写成 ask——笔误）
tri1 = T * mask                        # 第 3 行：边界 |x|=w 处 T=0，无半值歧义
print(f"[3.1] 与 base.tri(x/w) 最大差异 = "
      f"{np.max(np.abs(tri1 - tri(x / w))):.1f}（同一函数的两种写法）")

TRI = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(tri1))) * dx   # 居中 FFT × dx
TRI_an = w * sinc(w * fx) ** 2         # tri(x/w) <-> w sinc^2(w fx)

print(f"[3.1] 三角函数 FFT 峰值 = {np.abs(TRI).max():.6f}"
      f"（解析 w = 0.1 = 三角形面积）")
in_main = np.abs(fx) <= 30
err_tri = np.max(np.abs(np.abs(TRI[in_main]) - TRI_an[in_main]))
print(f"[3.1] FFT 与解析最大差异（|fx|<=30 内）= {err_tri:.2e}")

# ---------- 习题 3.3：圆函数 ----------
L2, M2 = 0.2, 200
dx2 = L2 / M2
x2, fx2 = coords(L2, M2)
y2 = x2
X, Y = np.meshgrid(x2, y2)             # 'xy'：行 <-> y
R = 0.015                              # 圆半径 (m) = 15 个样本

# 两种编码只差"圆周上的样本"：书风格（<=，边界取 1）vs base.circle（<）
circ_le = (np.hypot(X, Y) <= R).astype(float)   # 习题 3.3(a) 的书风格编码
circ_lt = circle(X, Y, R)                       # base 参考实现（严格小于）
n_bd = int(np.sum((np.hypot(X, Y) <= R).astype(int)) - np.sum((np.hypot(X, Y) < R).astype(int)))
print(f"[3.3] 半径 R = {R} m，两编码相差的圆周样本数 = {n_bd}"
      f"（对本例网格：{n_bd} 个）")

G_le = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(circ_le))) * dx2 ** 2
G_lt = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(circ_lt))) * dx2 ** 2
row = M2 // 2                          # fy=0 剖面所在行
J_an = jinc(fx2, R)                    # 半径 R 圆盘的解析频谱（Goodman 约定锚点）
print(f"[3.3] 圆盘 FFT 峰值 = {np.abs(G_le).max():.6e}"
      f"（解析面积 πR² = {np.pi * R ** 2:.6e}）")
in_main2 = np.abs(fx2) <= 30
err_le = np.max(np.abs(np.abs(G_le[row, :][in_main2]) - np.abs(J_an[in_main2])))
err_lt = np.max(np.abs(np.abs(G_lt[row, :][in_main2]) - np.abs(J_an[in_main2])))
print(f"[3.3] fy=0 剖面与 jinc 最大差异（|fx|<=30 内）："
      f"边界含圆周 {err_le:.2e}；不含 {err_lt:.2e}（相对峰值 ~1%，圆盘阶梯离散的正常水平）")

# ---------- 画图 ----------
fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.0))

axes[0, 0].plot(x, tri1, "o-", ms=4, mfc="none")
axes[0, 0].set_xlim(-0.3, 0.3); axes[0, 0].set_ylim(0, 1.2)
axes[0, 0].set_xlabel("x (m)")
axes[0, 0].set_title("Ex 3.1: triangle, half-base 0.1 m")

axes[0, 1].plot(fx, np.abs(TRI), label="discrete")
axes[0, 1].plot(fx, TRI_an, ":", label="analytic")
axes[0, 1].set_xlim(-60, 60); axes[0, 1].set_xlabel("fx (cyc/m)")
axes[0, 1].set_title(r"Ex 3.1: $|FFT|$ vs $w\,\mathrm{sinc}^2(w f_x)$")
axes[0, 1].legend(fontsize=8)

axes[1, 0].imshow(circ_le, extent=[x2[0], x2[-1], y2[0], y2[-1]],
                  origin="lower", cmap="gray")
axes[1, 0].set_xlabel("x (m)"); axes[1, 0].set_ylabel("y (m)")
axes[1, 0].set_title("Ex 3.3: circle, R = 0.015 m")

axes[1, 1].plot(fx2, np.abs(G_le[row, :]), "o", ms=4, mfc="none", label="discrete")
axes[1, 1].plot(fx2, np.abs(J_an), ":", label="analytic")
axes[1, 1].set_xlim(-60, 60); axes[1, 1].set_xlabel("fx (cyc/m)")
axes[1, 1].set_title(r"Ex 3.3: slice vs jinc $R J_1(2\pi R\rho)/\rho$")
axes[1, 1].legend(fontsize=8)

for ax in axes.flat:
    ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_exercises.png"), dpi=150)
print("已保存 out_exercises.png")
