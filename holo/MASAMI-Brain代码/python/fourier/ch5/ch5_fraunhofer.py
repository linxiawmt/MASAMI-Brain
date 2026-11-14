# -*- coding: utf-8 -*-
"""书 5.5 节 sqr_beam_FF 夫琅禾费传播例子(图 5.7 的 Python 复现)。

改动(相对 sqr_beam,书"Save As"清单):w = 0.011 m、用 propFF 传播、
观察面坐标按 L2 重算、nthroot(I2,3) 拉伸对比度(书上未列出的 I2 重算行
见笔记 warning)。NF = w^2/(lambda z) = 0.12 —— 夫琅禾费近似合理。
输出:out_fig5_7.png + 命令行数值检查(笔记 5.5 节"Python 实测")。

MATLAB/Python 差异:
  - propFF 返回 (u2, L2) 元组(MATLAB 是 [u2,L2] 输出列表);
  - MATLAB 的 nthroot(I2,3) 对应 np.cbrt(I2)(实数三次方根);
  - 书第 21 行注释 "obs ords" 应为 "obs coords"(排版小疵,见笔记 warning)。
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
from base import coords
from ch5_prop import propFF

# ---------- 书 sqr_beam 骨架(改动处见注释) ----------
L1 = 0.5                                             # 边长 (m)
M = 250                                              # 样本数
dx1 = L1 / M                                         # 源面采样间隔 = 2 mm
x1, _ = coords(L1, M)
y1 = x1

lam = 0.5e-6                                         # 波长 (m)
k = 2 * np.pi / lam                                  # 波数
w = 0.011                                            # 改动:源方孔半宽 (m)
z = 2000.0                                           # 传播距离 (m)

X1, Y1 = np.meshgrid(x1, y1)
u1 = np.where((np.abs(X1) < w) & (np.abs(Y1) < w), 1.0, 0.0)   # 式 (5.9)
I1 = np.abs(u1) ** 2                                 # 源面辐照度

# ---------- propFF 传播(观察面边长 L2 由式 (5.18) 决定) ----------
u2, L2 = propFF(u1, L1, lam, z)
dx2 = L2 / M                                         # 书改动清单:L2/M(= lam z/L1)
x2, _ = coords(L2, M)                                # 书改动清单:obs coords(原书误作 ords)
y2 = x2
I2 = np.abs(u2) ** 2                                 # 书省略的 I2 重算行(见笔记 warning)
NF = w ** 2 / (lam * z)

# ---------- 解析式 (5.22):I2 = (4w^2/lam z)^2 sinc^2(2w x2/lam z) sinc^2(...) ----------
I2_ana = (4 * w ** 2 / (lam * z)) ** 2 * \
    (np.sinc(2 * w * x2 / (lam * z)) ** 2)[:, None] * \
    (np.sinc(2 * w * y2 / (lam * z)) ** 2)[None, :]

# ---------- 图 5.7:(a) 图样 I2^(1/3) 拉伸;(b) 剖面 + 解析散点 ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
axs[0].imshow(np.cbrt(I2), extent=[x2[0], x2[-1] + dx2, y2[0], y2[-1] + dx2],
              origin="lower", cmap="gray")
axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)")
axs[0].set_title(f"z= {z:g} m"); axs[0].set_aspect("equal")
row = M // 2
axs[1].plot(x2, I2[row, :], label="propFF (250 pts)")
axs[1].plot(x2[::10], I2_ana[row, ::10], "o", mfc="none", label="analytic (5.22)")
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("Irradiance")
axs[1].set_title(f"z= {z:g} m"); axs[1].legend(); axs[1].grid(alpha=0.3)
fig.tight_layout()
outdir = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(outdir, "out_fig5_7.png"), dpi=110)
plt.close(fig)

# ---------- 数值检查(笔记 5.5 节"Python 实测") ----------
peak_ana = (4 * w ** 2 / (lam * z)) ** 2             # 解析峰值
zero_ana = lam * z / (2 * w)                         # sinc 第一零点位置
side_ana = 0.045 * peak_ana                          # 第一旁瓣(解析近似值)
i0 = np.argmin(np.abs(x2 - zero_ana))                # 最近零点的样本
side_num = np.max(I2[row, i0:i0 + int(round(zero_ana / dx2))])   # 第一、二零点之间
P1 = np.sum(I1) * dx1 ** 2                           # 源面功率 = (2w)^2
P2 = np.sum(I2) * dx2 ** 2                           # 观察面功率(定标各自正确)

print("== 书 5.5 节 sqr_beam_FF 数值检查 ==")
print(f"菲涅耳数 NF                 = {NF:.3f}   (书 0.12)")
print(f"L2 = lam z/dx1              = {L2:g} m , dx2 = {dx2:g} m (临界采样: L2 = L1)")
print(f"I2 峰值 解析 (4w^2/lam z)^2 = {peak_ana:.4f}")
print(f"I2 峰值 实测                = {I2[row, row]:.4f}")
print(f"sinc 第一零点 解析          = {zero_ana:.4f} m, 实测最近样本 = {x2[i0]:.4f} m")
print(f"第一旁瓣 解析 ~0.045*峰值   = {side_ana:.4f} , 实测 = {side_num:.4f}")
print(f"总功率 源面 = {P1:.3e} , 观察面 = {P2:.3e}")
print("已保存 out_fig5_7.png")
