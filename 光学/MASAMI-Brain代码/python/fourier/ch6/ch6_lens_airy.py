# -*- coding: utf-8 -*-
"""书 6.3 节 / 习题 6.4:f/10 理想透镜焦面的 Airy 图样。

复现:图 6.7(四次方根拉伸的焦面辐照度 + x 剖面)。
验证:
  (1) 焦面场 = propFF(z=f)(式 6.18)——直接用第 5 章的 propFF;
  (2) 解析 Airy 图样 I(r) = I0·[2J1(v)/v]²,v = k·wL·r/f,
      I0 = (π·wL²/(λf))²(单位振幅平面波全孔径聚焦);
  (3) 第一零点 1.22·λf/D(习题 6.4(a) 的解析推导);
  (4) 习题 6.4(b):缩小瞳直径,焦斑按比例变宽。
输出:out_fig6_7.png。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import j1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "ch5"))
from base import coords, circle                # noqa: E402  base.circle:半径 R 圆盘
from ch5_prop import propFF                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- 习题 6.4(a) 参数:f/10,2wL=25 mm ----------
lam = 0.5e-6
f = 0.25                                       # 焦距 250 mm
wL = 0.0125                                    # 孔径半径(2wL = 25 mm → f/10)
M = 250
L = 0.25                                       # 阵列边长 250 mm
dx = L / M
x, _ = coords(L, M)
X, Y = np.meshgrid(x, x)

# 式 (6.13):P = circ(√(x²+y²)/wL) —— Goodman 约定,半径 = wL。
# base.circle(X, Y, R) 画半径 R 圆盘,故 R = wL(等价 MATLAB: sqrt(X.^2+Y.^2)<=wL)。
u1 = circle(X, Y, wL)                          # 单位振幅平面波 × 瞳

u2, L2 = propFF(u1, L, lam, f)                 # 式 (6.18):z 换成 f
dx2 = L2 / M
x2, _ = coords(L2, M)
I2 = np.abs(u2) ** 2
prof = M // 2

# ---- 验证 (2):解析 Airy ----
lf = lam * f
I0 = (np.pi * wL ** 2 / lf) ** 2               # 理论峰值
v = np.pi * 2 * wL * x2 / lf                   # v = k wL r / f = πD r/(λf)
airy = I0 * np.where(np.abs(v) < 1e-9, 1.0,
                     (2 * j1(np.where(np.abs(v) < 1e-9, 1.0, v)) /
                      np.where(np.abs(v) < 1e-9, 1.0, v)) ** 2)
rel_prof = np.max(np.abs(I2[prof, :] - airy)) / I0
print(f"[2] 峰值:解析 I0 = {I0:.3e},数值 = {I2.max():.3e}"
      f"(差 {abs(I2.max()-I0)/I0*100:.3f}%);剖面最大相对差 = {rel_prof*100:.3f}%")

# ---- 验证 (3):第一零点 1.22 λf/D ----
# 注:0.5 um 网格采不到解析零点(最近样本幅值 ~2e-2),故找首个局部极小
zero_analytic = 1.22 * lf / (2 * wL)
Pk = np.abs(u2[prof, :])
ic = int(np.argmax(Pk))


def first_null(prof1d, ic):
    """从峰位向两侧走,找第一个严格局部极小(第一暗环的最近样本)。"""
    out = []
    for step in (+1, -1):
        i = ic
        while 0 < i + step < len(prof1d) - 1 and prof1d[i + step] < prof1d[i]:
            i += step
        out.append(x2[i])          # 停在最后一个"仍在下降"的样本 = 局部极小
    return out


nulls = first_null(Pk, ic)
print(f"[3] 第一零点:解析 1.22λf/D = {zero_analytic*1e6:.3f} um,"
      f"数值最近样本 = (+{abs(nulls[0])*1e6:.3f}, {nulls[1]*1e6:.3f}) um")

# ---- 验证 (4):习题 6.4(b) 缩小瞳(2wL = 12.5 mm)----
u1b = circle(X, Y, wL / 2)
I2b = np.abs(propFF(u1b, L, lam, f)[0]) ** 2
Pk_b = I2b[prof, :]
nulls_b = first_null(Pk_b, int(np.argmax(Pk_b)))
print(f"[4] 瞳减半:第一零点解析 {2*zero_analytic*1e6:.3f} um(翻倍),"
      f"数值最近样本 = (+{abs(nulls_b[0])*1e6:.3f}, {nulls_b[1]*1e6:.3f}) um;"
      f"峰值降为 {(Pk_b.max()/I2.max()):.4f}"
      f"(面积比 1/4,再摊到 2 倍宽的主瓣 → 峰高 1/16)")

# ---- 图 6.7 复现 ----
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
axs[0].imshow(I2 ** 0.25, origin="lower",
              extent=[x2[0], x2[-1], x2[0], x2[-1]], cmap="gray")
axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)")
axs[0].set_title("I2^(1/4) (fourth-root display)")
axs[1].plot(x2 * 1e5, I2[prof, :], label="discrete (propFF)")
axs[1].plot(x2 * 1e5, airy, "--", label="analytic Airy")
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("Irradiance")
axs[1].set_xlim(-2, 2); axs[1].legend()
axs[1].set_title("f/10 lens, focal plane (cf. Fig. 6.7)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_7.png"), dpi=150)
plt.close(fig)
print("done ch6_lens_airy")
