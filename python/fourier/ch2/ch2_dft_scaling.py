# ===== 书 2.3 / 2.5 节：DFT 定标与周期延拓（式 2.19-2.26，复现图 2.6） =====
# 书图 2.6 参数：L = 20, Δx = 1, N = 20 个样本，f(x) = rect(x/5)（半宽 2.5）。
# 验证三件事：
#   (1) 式 (2.19) 的 DFT 定义（含 exp[-j2π(pm/M+qn/N)] 核）与 np.fft.fft 逐元素一致
#       （前提：数据先排成"移位排列"——index 0 对应坐标 0，见 2.4 节）；
#   (2) 乘上 Δx 后，FFT 样本"跟随"周期延拓谱 F_P（Poisson 求和），
#       而与解析谱 F 的偏差 = 混叠（书 2.5 节图 2.6(b) 的高频瓣略大）；
#   (3) 式 (2.19)+(2.20) 连用（1/MN 因子）恢复原数组（傅里叶积分定理的离散版）。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import rect, coords

L = 20.0          # 采样空间边长 [m]（书图 2.6：L=20）
dx = 1.0          # 采样间隔 [m]
N = int(round(L / dx))
x, f = coords(L, N)          # 居中坐标：x ∈ [-10, 9]，f ∈ [-0.5, 0.4]
g = rect(x / 5.0)            # rect 半宽 2.5 → 样本 x=-2..2 共 5 个为 1，面积 5

# ---------- (1) 式 (2.19) 矩阵 DFT vs np.fft.fft（都在移位排列下） ----------
g_shift = np.fft.ifftshift(g)                     # 居中 → 移位（index 0 = 坐标 0）
m = np.arange(N)                                  # 式 (2.19) 的 m = 0..N-1（移位排列）
p = np.arange(N)[:, None]
W = np.exp(-2j * np.pi * p * m / N)               # DFT 核 exp[-j2π pm/M]
G_matrix = W @ g_shift                            # 式 (2.19)（无 Δx 因子）
G_fft = np.fft.fft(g_shift)
err_fft = np.max(np.abs(G_matrix - G_fft))
print(f"[1] 式(2.19)矩阵 DFT 与 np.fft.fft 最大差异 = {err_fft:.2e}")

# ---------- (2) FFT×Δx 样本 = 周期延拓谱 F_P 的样本 ----------
G = np.fft.fftshift(G_fft) * dx                   # 乘 Δx 回连续量纲，再居中排列
F_analytic = 5.0 * np.sinc(5.0 * f)               # 解析谱：rect(x/5) ↔ 5 sinc(5f)
k = np.arange(-3, 4)                              # 周期延拓：F_P(f) = Σ_k F(f - k/Δx)
F_periodic = np.sum([5.0 * np.sinc(5.0 * (f - kk / dx)) for kk in k], axis=0)
err_per = np.max(np.abs(G - F_periodic))
err_alias = np.max(np.abs(G - F_analytic))
print(f"[2] FFT×Δx 与周期谱 F_P 最大差异 = {err_per:.2e}")
print("    （FFT 样本只是'跟随'周期谱：空域矩形窗把 F_s 平滑了一层（式 2.24 的 sinc 卷积），")
print("      加上样本化效应，残差 ~1e-2；书图 2.6(b) 的点与虚线同样只差一点点）")
print(f"    FFT×Δx 与解析谱 F  最大差异 = {err_alias:.3e}（这就是混叠主项，瓣高处偏大）")
print(f"    DC 样本 G(0) = {G[N // 2].real:.6f}（= 中心纵坐标 ∫g dx = 5）")

# ---------- (3) 逆变换恢复原数组（式 2.20 的 1/MN 因子在 ifft 内置） ----------
g_back = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(G / dx)))
err_inv = np.max(np.abs(g_back - g))
print(f"[3] DFT → DFT^-1 往返恢复原数组最大差异 = {err_inv:.2e}")

# ---------- 画图（复现书图 2.6） ----------
fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2))
axes[0].plot(x, rect(x / 5.0), "k-", lw=1, label="analytic")
axes[0].plot(x - L, rect((x - L) / 5.0), "k--", lw=0.8, alpha=0.6)
axes[0].plot(x + L, rect((x + L) / 5.0), "k--", lw=0.8, alpha=0.6)
axes[0].plot(x, g, "o", ms=6, color="crimson", mfc="none", label="sampled (N=20)")
axes[0].set_xlim(-30, 30); axes[0].set_ylim(-0.05, 1.3)
axes[0].set_xlabel("x"); axes[0].set_ylabel("f(x)")
axes[0].set_title("Periodic extension: rect, L = 20, dx = 1")
axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.3)

axes[1].plot(f, np.abs(F_analytic), "k-", lw=1.2, label="analytic F")
axes[1].plot(f, np.abs(F_periodic), "k--", lw=0.9, alpha=0.7, label="periodic F_P")
axes[1].plot(f, np.abs(G), "o", ms=6, color="crimson", mfc="none", label="FFT x dx (N=20)")
axes[1].set_xlim(-1.5, 1.5); axes[1].set_ylim(0, 5.6)
axes[1].set_xlabel(r"$f_X$"); axes[1].set_ylabel(r"$|F(f_X)|$")
axes[1].set_title("Spectrum: FFT samples follow the periodic spectrum (aliasing lifts the sidelobes)")
axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_periodic_ext.png"), dpi=150)
print("已保存 out_periodic_ext.png")
