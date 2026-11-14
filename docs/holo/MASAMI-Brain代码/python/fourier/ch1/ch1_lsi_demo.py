# ===== 书 1.6 节：线性空间不变系统 —— 卷积定理数值演示（式 1.11-1.13） =====
# 线性空间不变系统的输出是卷积：  g2(x,y) = (g1 * h)(x,y)
# 两条等价路线：
#   A. 直接卷积（对应式 1.11 的积分，物理直观）
#   B. 频域法（式 1.13）：FFT(g1) * H -> IFFT，H 是传递函数
# 本脚本用"零填充 + FFT"实现线性卷积，与 convolve2d 的直接结果对比，
# 两者应只差浮点舍入（~1e-16）。
# 注意：不填充直接用 fft2 相乘得到的是"循环卷积"，边界会互相缠绕——
# 这是第 2 章"计算以外的世界"和第 5 章角谱法的核心话题，这里先埋个伏笔。
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d

N = 256            # 物面网格每边点数
L = 1.0            # 物面视场 [m]
dx = L / N
x = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(x, x)

# ---------- 输入物 g1：圆盘 + 小方块 + 双缝（简单"物体"） ----------
g1 = np.where(np.hypot(X + 0.25, Y + 0.22) < 0.10, 1.0, 0.0)      # 圆盘
g1 += np.where((np.abs(X - 0.22) < 0.05) & (np.abs(Y + 0.20) < 0.05), 1.0, 0.0)  # 方块
g1 += np.where((np.abs(X - 0.05) < 0.012) & (np.abs(Y - 0.25) < 0.10), 1.0, 0.0)  # 细棒 1
g1 += np.where((np.abs(X + 0.05) < 0.012) & (np.abs(Y - 0.25) < 0.10), 1.0, 0.0)  # 细棒 2
g1 = np.clip(g1, 0, 1)

# ---------- 冲激响应 h：小圆盘（"模糊斑"，即光学里的点扩散函数 PSF） ----------
rh = 0.03
Kh = 31                       # h 网格边长（奇数，中心对齐）
xh = (np.arange(Kh) - Kh // 2) * dx
Xh, Yh = np.meshgrid(xh, xh)
h = np.where(np.hypot(Xh, Yh) < rh, 1.0, 0.0)
h /= h.sum()                  # 归一化能量：直流增益 H(0,0)=1，像不整体变亮

# ---------- 路线 B：频域法（零填充到 N+Kh-1，做线性卷积） ----------
P = N + Kh - 1
G2freq = np.fft.ifft2(np.fft.fft2(g1, s=(P, P)) * np.fft.fft2(h, s=(P, P)))
g2_freq = G2freq[(Kh // 2):(Kh // 2) + N, (Kh // 2):(Kh // 2) + N].real
# 裁切说明：h 中心在 (0,0)，其样本平移了 Kh//2；线性卷积结果里
# g1 的第 u 个样本对应整体输出的第 u + Kh//2 个。

# ---------- 路线 A：直接二维卷积（scipy 实现式 1.11 的求和） ----------
g2_direct = convolve2d(g1, h, mode="same", boundary="fill")

diff = np.max(np.abs(g2_direct - g2_freq))
rel = np.abs(g2_direct - g2_freq) / g2_freq.max()   # 相对误差（显示用）
print(f"两种路线结果最大差异 = {diff:.3e}（相对量 {rel.max():.1e}，浮点舍入量级）")
print(f"g2 峰值 = {g2_freq.max():.4f}（模糊后仍 > 单点峰值，因为物有连通结构）")

# ---------- 画图 ----------
fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
extent = [x[0], x[-1], x[0], x[-1]]
h_full = np.zeros_like(g1)               # 把 h 镶进全画面，统一显示比例
o = (N - Kh) // 2
h_full[o:o + Kh, o:o + Kh] = h
axes[0].imshow(g1, cmap="gray", extent=extent, origin="lower")
axes[0].set_title("input  g1 (object)")
axes[1].imshow(h_full, cmap="gray", extent=extent, origin="lower")
axes[1].set_title(f"impulse response h (PSF), {Kh}x{Kh} grid")
axes[2].imshow(g2_freq, cmap="gray", extent=extent, origin="lower")
axes[2].set_title("output  g2 = g1 * h")
im = axes[3].imshow(rel, extent=extent, origin="lower")
axes[3].set_title(f"|direct - spectral| / max(g2)\n(max {rel.max():.0e})")
plt.colorbar(im, ax=axes[3], fraction=0.046)
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Linear space-invariant system: convolution vs. transfer-function route")
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(__file__.replace("ch1_lsi_demo.py", "out_lsi_demo.png"), dpi=150)
print("已保存 out_lsi_demo.png")
