# -*- coding: utf-8 -*-
"""书 6.2 节:sqr_beam_focus 聚焦方形光束 + 式 (6.10)/(6.11) 数值验证。

复现:图 6.5(zf = z = 2000 m 的聚焦光斑与剖面)。
验证:
  (1) "聚焦图样 = 缩放夫琅禾费图样"(书 6.2 节原话,式 6.18 的预言):
      focus(zf=z) 后 propTF(z) 的辐照度 == propFF(u1, z=f) 的辐照度;
  (2) 式 (6.11) 边界:|zf| >= 671 m(B1 取书用的 98 cyc/m);
      zf=−671(发散,推过边界)/ +671(先聚焦,图样变小)对照。
输出:out_fig6_5.png、out_fig6_focus_bounds.png。
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
from base import coords                        # noqa: E402
from ch5_prop import propTF, propFF            # noqa: E402
from ch6_transmittance import focus            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- 书 5.3 节 sqr_beam 源面 ----------
L1 = 0.5
M = 250
dx1 = L1 / M
x1, _ = coords(L1, M)
y1 = x1
lam = 0.5e-6
w = 0.051
z = 2000.0
X1, Y1 = np.meshgrid(x1, y1)
u1 = np.where((np.abs(X1) < w) & (np.abs(Y1) < w), 1.0, 0.0)

# ---------- 书 6.2 节插入块 ----------
zf = 2000.0
u1f = focus(u1, L1, lam, zf)
u2 = propTF(u1f, L1, lam, z)                   # 传播
I2 = np.abs(u2) ** 2
prof = M // 2                                  # y=0 行

# ---- 验证 (1):与 propFF(z=f) 对拍(式 6.18:焦面场 = 夫琅禾费图样) ----
u2ff, L2 = propFF(u1, L1, lam, zf)             # propFF 的观察面 L2 = λz/Δx1 = 0.5 m
I2ff = np.abs(u2ff) ** 2
rel = np.max(np.abs(I2 - I2ff)) / I2.max()
print(f"[1] focus(zf=z)+propTF vs propFF(z=f):辐照度最大相对差 = {rel:.2e}"
      f"(L2 = {L2} m,与源面同网格)")

# ---- 图 6.5 复现 ----
fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.8))
axs[0].imshow(I2, origin="lower", extent=[x1[0], x1[-1], y1[0], y1[-1]], cmap="gray")
axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)"); axs[0].set_title("z= 2000 m")
axs[1].plot(x1, I2[prof, :])
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("Irradiance"); axs[1].set_title("z= 2000 m")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_5.png"), dpi=150); plt.close(fig)
prof1d = I2[prof, :]
above = x1[prof1d > prof1d.max() / 2]
print(f"[1b] 焦斑峰值 I2max = {I2.max():.1f},半高宽 ≈ {(above[-1]-above[0])*1e3:.2f} mm(书图 ~mm 量级)")

# ---- 验证 (2):式 (6.11) 边界 |zf| >= 671 m ----
B1_book = 5 / w                                # 书用 98(实为 5/w;式 5.10 原值 117.3)
D1 = 2 * w
zf_bound = (D1 / 2 / lam) / (1 / (2 * dx1) - B1_book)     # 式 (6.11)
print(f"[2] B1(书用)= {B1_book:.1f};式(6.11) 边界 |zf| = {zf_bound:.1f} m(书 671 m)")
print(f"    若用 B1=117.3:|zf| = {(D1/2/lam)/(1/(2*dx1)-5.98/w):.1f} m")

rows = []
for zft in (+zf_bound, -zf_bound):
    u2b = propTF(focus(u1, L1, lam, zft), L1, lam, z)
    rows.append(np.abs(u2b)[prof, :])
    print(f"    zf={zft:+.0f}:|u2| 峰 {np.abs(u2b).max():.3f},边缘均值 "
          f"{0.5*(np.abs(u2b)[prof,:5].mean()+np.abs(u2b)[prof,-5:].mean()):.3f}")

fig, axs = plt.subplots(2, 1, figsize=(7.2, 5.2), sharex=True)
axs[0].plot(x1, rows[0], color="seagreen")
axs[0].set_title(rf"zf=+{zf_bound:.0f} m (bound): smaller pattern, edge phase aliasing")
axs[1].plot(x1, rows[1], color="steelblue")
axs[1].set_title(rf"zf=$-${zf_bound:.0f} m (bound): diverging beam pushes over array")
for a in axs:
    a.set_ylabel("|u2|")
axs[1].set_xlabel("x (m)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_focus_bounds.png"),
                                dpi=150); plt.close(fig)
print("done ch6_sqr_beam_focus")
