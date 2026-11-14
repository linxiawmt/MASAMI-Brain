# ===== 书 2.6 节：周期卷积（式 2.27，复现图 2.7） =====
# 书图 2.7 参数：L = 20, Δx = 1, N = 20 个样本；
#   f = rect，支撑 D1 = 9（样本 x = -4..4 为 1）；h = 右三角形（上升斜坡），
#   支撑 D2 = 15（x ∈ [-7.5, 7.5]，h(x) = (x+7.5)/15）。
# D1 + D2 = 24 > L = 20 → 违反式 (2.27)：FFT 卷积 ≠ 解析卷积（边缘卷绕污染）。
# 演示两条修复路线：(a) 换更大的 L = 40 重新采样；(b) 同一组样本零填充后再 FFT。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import rect, coords


def ramp(x, D):
    """支撑宽度 D 的右三角形（上升斜坡，峰值 1 在右边缘）。"""
    return np.where(np.abs(x) <= D / 2, (x + D / 2) / D, 0.0)


def fft_conv(g, h, dx):
    """FFT 卷积（循环卷积，不填充）：输入居中排列，输出居中排列。"""
    G = np.fft.fft(np.fft.ifftshift(g)) * np.fft.fft(np.fft.ifftshift(h)) * dx
    return np.fft.fftshift(np.fft.ifft(G))


# ---------- 解析参照：细网格直接卷积，再对齐到 xa 坐标 ----------
xa = np.arange(-30.0, 30.0 + 1e-9, 0.005)       # 中心 index 恰为 len//2 ↔ x=0
dd = xa[1] - xa[0]
conv_full = np.convolve(rect(xa / 9.0), ramp(xa, 15.0)) * dd   # 支撑 [-24, 24]
c0 = len(conv_full) // 2                         # conv_full[c0] ↔ t = 0
conv_ref = conv_full[c0 - len(xa) // 2: c0 + len(xa) // 2 + 1]  # 对齐到 xa

# ---------- L = 20：违反式 (2.27) ----------
L, N = 20.0, 20
x, _ = coords(L, N)
g = rect(x / 9.0)            # D1 = 9
h = ramp(x, 15.0)            # D2 = 15
conv_fft = fft_conv(g, h, dx=1.0).real

err20 = np.max(np.abs(np.interp(x, xa, conv_ref) - conv_fft))
pk = np.argmax(conv_ref)
print(f"L=20（D1+D2=24 > L=20，违反式 2.27）：FFT 卷积与解析最大差异 = {err20:.3f}")
print(f"  解析峰值 = {conv_ref[pk]:.3f} @ x = {xa[pk]:.2f}；FFT 峰值 = {conv_fft.max():.3f}")
print("  （书图 2.7(d)：FFT 样本在高频端/边缘跟着周期卷积走，卷绕区多出一块错误面积）")

# ---------- 修复 1：L = 40 重新采样（D1+D2=24 < 40，满足式 2.27） ----------
L2, N2 = 40.0, 40
x2, _ = coords(L2, N2)
conv_fft2 = fft_conv(rect(x2 / 9.0), ramp(x2, 15.0), dx=1.0).real
err40 = np.max(np.abs(np.interp(x2, xa, conv_ref) - conv_fft2))
print(f"L=40（满足式 2.27）：FFT 卷积与解析最大差异 = {err40:.3e}（矩形求和近似误差量级）")

# ---------- 修复 2：零填充（同一组 N=20 样本，填到 N1+N2-1 再 FFT） ----------
P = N + N - 1                                 # 20+20-1 = 39（奇数：中心样本唯一）
g_pad = np.zeros(P); h_pad = np.zeros(P)
o = (P - N) // 2
g_pad[o:o + N] = g; h_pad[o:o + N] = h
Gp = np.fft.fft(np.fft.ifftshift(g_pad)) * np.fft.fft(np.fft.ifftshift(h_pad)) * 1.0
conv_pad = np.fft.fftshift(np.fft.ifft(Gp)).real
xp = (np.arange(P) - P // 2) * 1.0            # 填充网格的居中坐标
err_pad = np.max(np.abs(np.interp(xp, xa, conv_ref) - conv_pad))
print(f"零填充 N1+N2-1={P}：FFT 卷积与解析最大差异 = {err_pad:.3e}")

# ---------- 画图：复现图 2.7(d) + 两种修复 ----------
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.0))
axes[0].plot(xa, conv_ref, "k-", lw=1.4, label="analytic")
axes[0].plot(x, conv_fft, "o", ms=6, color="crimson", mfc="none", label="FFT (L=20)")
axes[0].set_title(f"L = 20: D1+D2 = 24 > L\nperiodic wrap-around, max err {err20:.2f}")
axes[1].plot(xa, conv_ref, "k-", lw=1.4, label="analytic")
axes[1].plot(x2, conv_fft2, "o", ms=5, color="teal", mfc="none", label="FFT (L=40)")
axes[1].set_title(f"L = 40: D1+D2 < L satisfied\nmax err {err40:.1e}")
axes[2].plot(xa, conv_ref, "k-", lw=1.4, label="analytic")
axes[2].plot(xp, conv_pad, "o", ms=5, color="darkorange", mfc="none",
             label="zero-padded FFT")
axes[2].set_title(f"zero-padded to N1+N2-1 = {P}\nmax err {err_pad:.1e}")
for ax in axes:
    ax.set_xlim(-30, 30); ax.set_ylim(-0.3, 7); ax.grid(True, alpha=0.3)
    ax.set_xlabel("x"); ax.set_ylabel("(f*h)(x)"); ax.legend(fontsize=8, loc="upper left")
fig.suptitle("Periodic convolution (book Fig. 2.7): rect D1=9 * right-triangle D2=15")
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_periodic_conv.png"), dpi=150)
print("已保存 out_periodic_conv.png")
