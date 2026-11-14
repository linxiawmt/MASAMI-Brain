# -*- coding: utf-8 -*-
"""书 6.1 节:sqr_beam_tilt 倾斜方形光束 + 采样判据数值验证。

复现:图 6.3(α=5e-5 rad、θ=45°的观察面辐照度)。
验证:
  (1) 频谱平移:tilt 后源谱峰移到 (cosθ,sinθ)·tanα/λ(式 6.4/6.5);
  (2) 观察面质心平移 r = z·tanα(图 6.2);
  (3) 判据试验(书 6.1 节末,θ=0):
      α=7.6e-5(书用 B1=98 得到的式 6.6 边界)→ 光束贴近边缘 + wrap-around;
      α=λ/2Δx1=1.25e-4 → 图案被网格边界切成两半。
输出:out_fig6_3.png、out_fig6_tilt_lim.png(保存在本目录)。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "ch5"))
from base import coords, ft2_centered          # noqa: E402  公共积木(ch1/ch2 起)
from ch5_prop import propTF                    # noqa: E402  第 5 章传播器
from ch6_transmittance import tilt             # noqa: E402  本章透过率(ch6_transmittance.py)

# ---------- 书 5.3 节 sqr_beam 的源面(与第 5 章相同) ----------
L1 = 0.5                                       # 边长 (m)
M = 250                                        # 样本数
dx1 = L1 / M                                   # 采样间隔 2 mm
x1, _ = coords(L1, M)
y1 = x1

lam = 0.5e-6                                   # 波长 (m)
w = 0.051                                      # 源方孔半宽 (m)
z = 2000.0                                     # 传播距离 (m)

X1, Y1 = np.meshgrid(x1, y1)
u1 = np.where((np.abs(X1) < w) & (np.abs(Y1) < w), 1.0, 0.0)   # rect 掩码

# ---------- 书 6.1 节插入块:alpha=5.0e-5, theta=45 deg ----------
deg = np.pi / 180
alpha = 5.0e-5                                 # rad
theta = 45 * deg
u1t = tilt(u1, L1, lam, alpha, theta)

# ---- 验证 (1):频谱平移量 = tanα/λ = 100 cyc/m,沿 45°(式 6.5) ----
G1 = ft2_centered(u1, dx1)                     # 未倾斜谱(居中)
G1t = ft2_centered(u1t, dx1)                   # 倾斜谱(居中)
peak = np.unravel_index(np.argmax(np.abs(G1t)), G1t.shape)
f_axis = np.fft.fftshift(np.fft.fftfreq(M, dx1))   # 与 base.coords 的 f 轴一致
f0 = np.tan(alpha) / lam                       # 理论平移 100 cyc/m
fx_pk, fy_pk = f_axis[peak[1]], f_axis[peak[0]]
print(f"[1] 谱峰位置 = ({fx_pk:.2f}, {fy_pk:.2f}) cyc/m, "
      f"理论 = ({f0*np.cos(theta):.2f}, {f0*np.sin(theta):.2f}) cyc/m")

# ---- 验证 (2):观察面质心平移 r = z tanα = 0.1 m,45° 方向 ----
u2 = propTF(u1t, L1, lam, z)                   # 书 6.1 节:propTF
I2 = np.abs(u2) ** 2
P2 = I2.sum()
xc = (I2.sum(axis=0) * x1).sum() / P2          # 质心 x(能量加权)
yc = (I2.sum(axis=1) * y1).sum() / P2
print(f"[2] 观察面质心 = ({xc*1e3:.2f}, {yc*1e3:.2f}) mm, "
      f"理论 = ({z*np.tan(alpha)*np.cos(theta)*1e3:.2f}, "
      f"{z*np.tan(alpha)*np.sin(theta)*1e3:.2f}) mm")

# ---- 图 6.3 复现 ----
fig, ax = plt.subplots(figsize=(4.6, 4.0))
ax.imshow(I2, origin="lower", extent=[x1[0], x1[-1], y1[0], y1[-1]],
          cmap="gray", vmax=I2.max() * 0.4)
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_title("z= 2000 m")
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "out_fig6_3.png"), dpi=150); plt.close(fig)

# ---- 验证 (3):判据边界试验(θ=0,书 6.1 节末) ----
B1_book = 5 / w                                # 书 6.1 节用的 B1 ≈ 98(实为 5/w,见笔记勘误)
B1_598 = 5.98 / w                              # 式 (5.10) 原值 = 117.3(前后不一致)
a_crit = lam * (1 / (2 * dx1) - B1_book)       # 式 (6.6):7.6e-5 rad
print(f"[3] B1(书用)= {B1_book:.1f} cyc/m(5/w);式(5.10)原值 {B1_598:.1f};"
      f"式(6.6) 边界 alpha = {a_crit:.2e} rad(书 7.6e-5)")
print(f"    若用 B1=117.3:边界 = {lam*(1/(2*dx1)-B1_598):.2e} rad")

u_edge = tilt(u1, L1, lam, a_crit, 0.0)        # 贴边:能量到阵列边缘
u2e = propTF(u_edge, L1, lam, z)
u_split = tilt(u1, L1, lam, lam / (2 * dx1), 0.0)   # α=λ/2Δx1=1.25e-4
u2s = propTF(u_split, L1, lam, z)
prof = M // 2                                  # y=0 行(numpy 0-based)

fig, axs = plt.subplots(3, 1, figsize=(7.2, 6.4), sharex=True)
axs[0].plot(x1, np.abs(propTF(tilt(u1, L1, lam, 0.0, 0.0), L1, lam, z))[prof, :])
axs[0].set_title(r"no tilt ($\alpha$=0)")
axs[1].plot(x1, np.abs(u2e)[prof, :], color="crimson")
axs[1].set_title(rf"$\alpha$={a_crit:.2e} (Eq. 6.6 bound, B1=98): beam at edge + wrap-around")
axs[2].plot(x1, np.abs(u2s)[prof, :], color="darkorange")
axs[2].set_title(rf"$\alpha$=$\lambda/2\Delta x_1$={lam/(2*dx1):.2e}: pattern splits across boundary")
for a in axs:
    a.set_ylabel("|u2|")
axs[2].set_xlabel("x (m)")
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "out_fig6_tilt_lim.png"), dpi=150); plt.close(fig)

# wrap-around 定量:贴边情形对侧(负 x)边缘的能量
neg_edge = np.abs(u2e)[prof, :M // 4].sum()
pos_edge = np.abs(u2e)[prof, 3 * M // 4:].sum()
print(f"[4] alpha=边界值:对侧(-x)边缘能量占比 {neg_edge/np.abs(u2e)[prof,:].sum():.3f},"
      f"本侧(+x)边缘 {pos_edge/np.abs(u2e)[prof,:].sum():.3f}(wrap-around 证据)")
print("done ch6_sqr_beam_tilt")
