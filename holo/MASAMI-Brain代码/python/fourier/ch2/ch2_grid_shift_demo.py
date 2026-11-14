# ===== 书 2.4 节：坐标、索引、居中与移位（式 2.21/2.22，复现图 2.4 / 2.5） =====
# 三种排列（书图 2.4）：
#   analytic ：指标带正负（m = -M/2 .. M/2-1），纸面公式用；
#   centered ：放进计算机后的排列，零坐标在 index M/2+1（MATLAB 1-based）
#              ↔ Python 0-based 的 index N//2；
#   shifted  ：FFT 约定 index 1（MATLAB）/ 0（Python）= 零坐标，FFT 前必须转成它。
# centered ↔ shifted 用 fftshift / ifftshift 互换（偶数 N 两者相同，奇数 N 不同！）。
# 本脚本用 N=8 小网格把三种排列打印出来，并复现图 2.4（1D）与图 2.5（2D 象限交换）。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import rect, coords

# ---------- 1D：N=8 小网格，rect 半宽 1.5（3 个样本为 1） ----------
N, L = 8, 8.0
dx = L / N
x, f = coords(L, N)                        # 式 (2.21)/(2.22)：居中升序坐标
g_centered = rect(x / 3.0)                 # 居中排列：|x|<1.5 → x=-1,0,1
g_shifted = np.fft.ifftshift(g_centered)   # 居中 → 移位（FFT 约定）

print("N=8, L=8, dx=1 的三种排列（数值对应书图 2.4）")
print("MATLAB 1-based index:", np.arange(1, N + 1))
print("Python 0-based index:", np.arange(N))
print("coords x (式 2.21)   :", x)
print("g 居中 centered      :", g_centered.astype(int))
print("g 移位 shifted(FFT用):", g_shifted.astype(int))
print("freq 坐标 f (式 2.22):", f)
print("roundtrip fftshift(ifftshift(g)) == g :",
      np.allclose(np.fft.fftshift(g_shifted), g_centered))
print("注意：偶数 N 时 fftshift 与 ifftshift 等价；奇数 N 不等价（见笔记坑与陷阱）")

# ---------- 2D：象限交换（书图 2.5） ----------
N2 = 8
x2, _ = coords(L, N2)
X, Y = np.meshgrid(x2, x2)                # 'xy'：行 ↔ y，列 ↔ x（与 MATLAB 一致）
A_centered = (np.hypot(X, Y) < 1.5).astype(float)   # 居中圆盘，[0,0] 在 (N/2, M/2)
A_shifted = np.fft.ifftshift(A_centered)            # [0,0] 挪到 (0, 0)
print("\n2D 居中排列零坐标位置 (Python 行,列) =",
      tuple(np.argwhere(A_centered > 0).mean(axis=0).round(1)),
      "；MATLAB (i,j) 记法 = (N/2+1, M/2+1) = (5, 5)")
print("2D 移位排列零坐标位置 (Python 行,列) =", tuple(np.argwhere(A_shifted > 0)[0]))

# ---------- 画图：上排复现图 2.4 三种 1D 排列，下排复现图 2.5 象限交换 ----------
fig = plt.figure(figsize=(11.5, 7.6))
gs = fig.add_gridspec(2, 3)

ax = fig.add_subplot(gs[0, 0])
ax.plot(x, g_centered, "o-", ms=5, color="crimson", mfc="none")
ax.set_title("centered: zero coord at index N//2\n(apply ifftshift before fft)")
ax.set_xticks(x); ax.grid(True, alpha=0.3); ax.set_ylim(-0.2, 1.25)

ax = fig.add_subplot(gs[0, 1])
ax.plot(np.arange(N), g_shifted, "s-", ms=5, color="teal", mfc="none")
ax.set_title("shifted for FFT: zero coord at index 0")
ax.set_xticks(np.arange(N)); ax.grid(True, alpha=0.3); ax.set_ylim(-0.2, 1.25)

ax = fig.add_subplot(gs[0, 2])
G = np.fft.fft(g_shifted) * dx
ax.plot(f, np.abs(np.fft.fftshift(G)), "o-", ms=5, color="slateblue", mfc="none")
ax.set_title("|FFT x dx| back on centered freq axis")
ax.grid(True, alpha=0.3)

# 象限编号（书图 2.5：居中 = 3,4/2,1；移位 = 1,2/4,3）
q = np.zeros((N2, N2), int)
q[:N2 // 2, :N2 // 2] = 3; q[:N2 // 2, N2 // 2:] = 4
q[N2 // 2:, :N2 // 2] = 2; q[N2 // 2:, N2 // 2:] = 1
q_shift = np.fft.ifftshift(q)
for k, (A, ttl) in enumerate([(q, "2D centered (quadrant No.)\n[0,0] at MATLAB (N/2+1, M/2+1)"),
                              (q_shift, "2D shifted for FFT (quadrant No.)\n[0,0] at index (1,1)/(0,0)")]):
    ax = fig.add_subplot(gs[1, k])
    ax.imshow(A, cmap="Pastel1", origin="upper")
    for i in range(N2):
        for j in range(N2):
            ax.text(j, i, str(A[i, j]), ha="center", va="center", fontsize=9)
    ax.set_title(ttl); ax.tick_params(labelbottom=False, labelleft=False)

ax = fig.add_subplot(gs[1, 2])
ax.plot(np.abs(G), "s-", ms=5, color="darkorange", mfc="none")
ax.axvline(N2 // 2 - 0.5, color="gray", ls=":", lw=1)
ax.text(N2 // 2 - 0.2, 4.6, "raw FFT output:\nleft half = zero & positive freq,\nright half = negative freq", fontsize=8)
ax.set_title("why fftshift exists: raw fft\nordering vs. centered axis")
ax.set_xticks(np.arange(N2)); ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_arrangements.png"), dpi=150)
print("已保存 out_arrangements.png")
