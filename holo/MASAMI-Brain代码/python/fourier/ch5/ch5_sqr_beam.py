# -*- coding: utf-8 -*-
"""书 5.3 节 sqr_beam 方形光束传播例子(图 5.1-5.3 的 Python 复现)。

参数:L1 = 0.5 m、M = 250、dx1 = 2 mm、w = 0.051 m(跨孔径 51 样本)、
lambda = 0.5 um、z = 2000 m、NF = w^2/(lambda z) = 2.6(书:"reasonable")。
输出:out_fig5_1/2/3.png + 命令行数值检查(笔记 5.3 节"Python 实测数值检查")。

MATLAB/Python 差异:
  - MATLAB 偶数阵列的"零坐标"样本在 (M/2+1, M/2+1),numpy 在 (M//2, M//2),
    所以书第 35 行 I2(M/2+1,:) 对应 I2[M//2, :](M=250 时即第 126 行)。
  - imagesc + axis xy 对应 plt.imshow(origin='lower')。
  - MATLAB 基础版没有 rect(信号处理工具箱才有),书第 14 行用逻辑掩码等价实现。
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
from base import coords                              # 公共积木(ch1 起)
from ch5_prop import propTF, propIR                  # 本章三个传播器(ch5_prop.py)

# ---------- 书 sqr_beam 第 1-7 行:源面网格 ----------
L1 = 0.5                                             # 边长 (m)
M = 250                                              # 样本数
dx1 = L1 / M                                         # 源面采样间隔 = 2 mm
x1, _ = coords(L1, M)                                # 书第 6 行:-L1/2 : dx1 : L1/2-dx1
y1 = x1                                              # 书第 7 行

# ---------- 书第 8-14 行:方孔源场 ----------
lam = 0.5e-6                                         # 波长 (m)
k = 2 * np.pi / lam                                  # 波数
w = 0.051                                            # 源方孔半宽 (m)
z = 2000.0                                           # 传播距离 (m)

X1, Y1 = np.meshgrid(x1, y1)                         # 书第 13 行
u1 = np.where((np.abs(X1) < w) & (np.abs(Y1) < w), 1.0, 0.0)   # 书第 14 行
#     rect(X1/(2w))*rect(Y1/(2w));矩形的边界值歧义对均匀网格无影响
I1 = np.abs(u1) ** 2                                 # 书第 15 行:源面辐照度

# ---------- 书第 22 行:传播(TF 与 IR 各一次) ----------
u2_tf = propTF(u1, L1, lam, z)
u2_ir = propIR(u1, L1, lam, z)
I2_tf = np.abs(u2_tf) ** 2                           # 书第 26 行:观察面辐照度
row = M // 2                                         # y = 0 的行(numpy 0-based)
prof_tf = I2_tf[row, :]                              # 书第 35 行:I2(M/2+1, :)
x2 = x1                                              # 书第 24 行:观察面坐标 = 源面坐标
y2 = y1                                              # 书第 25 行

# ================= 以下为书外部分:出图 + 数值检查 =================
outdir = os.path.dirname(os.path.abspath(__file__))

# ---------- 图 5.1:源面辐照度(书 figure 1) ----------
fig, ax = plt.subplots(figsize=(4.6, 4.0))
ax.imshow(I1, extent=[x1[0], x1[-1] + dx1, y1[0], y1[-1] + dx1],
          origin="lower", cmap="gray")
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_title("z= 0 m")
ax.set_aspect("equal")
fig.tight_layout(); fig.savefig(os.path.join(outdir, "out_fig5_1.png"), dpi=110)
plt.close(fig)

# ---------- 图 5.2:观察面辐照度图样 + 剖面(书 figure 2/3) ----------
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
axs[0].imshow(I2_tf, extent=[x2[0], x2[-1] + dx1, y2[0], y2[-1] + dx1],
              origin="lower", cmap="gray")
axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)")
axs[0].set_title(f"z= {z:g} m"); axs[0].set_aspect("equal")
axs[1].plot(x2, prof_tf)
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("Irradiance")
axs[1].set_title(f"z= {z:g} m"); axs[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(outdir, "out_fig5_2.png"), dpi=110)
plt.close(fig)

# ---------- 图 5.3:幅度 + unwrap 相位剖面(书 figure 4/5) ----------
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
axs[0].plot(x2, np.abs(u2_tf[row, :]))
axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("Magnitude")
axs[0].set_title(f"z= {z:g} m"); axs[0].grid(alpha=0.3)
axs[1].plot(x2, np.unwrap(np.angle(u2_tf[row, :])))
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("Phase (rad)")
axs[1].set_title(f"z= {z:g} m"); axs[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(outdir, "out_fig5_3.png"), dpi=110)
plt.close(fig)

# ---------- 数值检查(笔记 5.3 节"Python 实测数值检查") ----------
I2_ir = np.abs(u2_ir) ** 2                           # IR 传播的观察面辐照度
n_ap = int(round(2 * w / dx1))                       # 跨孔径样本数
NF = w ** 2 / (lam * z)                              # 菲涅耳数
rel = np.max(np.abs(I2_tf - I2_ir)) / np.max(I2_tf)  # TF vs IR 最大相对差

c_idx = row                                          # x = 0 的样本(125)
p_side = I2_tf[row, c_idx + int(round(0.02 / dx1))]  # x = +0.02 m 侧峰
v_mid = I2_tf[row, c_idx + int(round(0.01 / dx1))]   # x = +0.01 m 谷
P1 = np.sum(I1) * dx1 ** 2                           # 式 (5.25):源面功率
P2 = np.sum(I2_tf) * dx1 ** 2                        # 观察面功率(TF)

print("== 书 5.3 节 sqr_beam 数值检查 ==")
print(f"跨孔径样本数 2w/dx1      = {n_ap}")
print(f"菲涅耳数 NF = w^2/(lam z) = {NF:.2f}")
print(f"TF 与 IR 辐照度最大相对差 = {rel:.1e}")
print(f"中央峰 I2(0,0)            = {I2_tf[row, c_idx]:.3f}")
print(f"x=+0.02 m 侧峰            = {p_side:.3f}")
print(f"x=+0.01 m 谷              = {v_mid:.3f}")
print(f"总功率 源面 = {P1:.6f} , 观察面 = {P2:.6f}")
print("已保存 out_fig5_1.png / out_fig5_2.png / out_fig5_3.png")
