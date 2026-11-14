# -*- coding: utf-8 -*-
"""书 5.7 节习题的 Python 数值验证:5.1 / 5.2 / 5.3 / 5.5 / 5.6 / 5.8。

5.4(RS 传播器)、5.7(三孔径夫琅禾费)、5.9(两步法)笔记只给提示,不在本脚本内。
输出:out_ex5_1.png(习题 5.1d 圆孔三距离)、out_ex5_5.png(习题 5.5a 高斯束)
+ 各题命令行数值(笔记习题节的"参考答案要点"数字全部来自本脚本实跑)。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfinv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
from base import coords
from ch5_prop import propTF, propIR

outdir = os.path.dirname(os.path.abspath(__file__))

# ======================================================================
# 习题 5.1 圆孔 + 临界采样:w = 0.05 m、lambda = 0.5 um、z = 1000 m、M = 500
# ======================================================================
print("=" * 60)
print("习题 5.1 圆孔临界采样")
w = 0.05                                             # 圆孔半径 (m)
lam = 0.5e-6                                         # 波长 (m)
M = 500                                              # 阵列 500x500
z0 = 1000.0                                          # 传播距离 (m)

L1 = np.sqrt(M * lam * z0)                           # 临界采样: L1^2 = M lam z(式 5.14)
dx = L1 / M                                          # 采样间隔
FN = 1 / (2 * dx)                                    # Nyquist 频率
B1 = 5.98 / w                                        # 源有效带宽(式 5.10 口径)
NF = w ** 2 / (lam * z0)                             # 菲涅耳数
print(f"(a) L1 = sqrt(M lam z) = {L1:g} m , dx = {dx:g} m , FN = {FN:g} cyc/m")
print(f"(b) B1 = 5.98/w = {B1:.1f} cyc/m < FN = {FN:g}? {B1 < FN} ; "
      f"跨直径样本 = {int(round(2*w/dx))}")
print(f"(c) NF = w^2/(lam z) = {NF:g}(NF=1 对应 z = {w**2/lam:g} m)")

x1, _ = coords(L1, M)
X1, Y1 = np.meshgrid(x1, x1)
u1 = np.where(X1 ** 2 + Y1 ** 2 < w ** 2, 1.0, 0.0)  # 离散圆盘(半径 w)
P_disk = np.sum(np.abs(u1) ** 2) * dx ** 2           # 离散圆盘功率
print(f"    离散圆盘功率 = {P_disk:.6f} vs 解析 pi w^2 = {np.pi*w**2:.6f}")

zs = [500.0, 1000.0, 2000.0]
fig, axs = plt.subplots(1, 3, figsize=(15.0, 5.0))
tag = ["TF over (2x)", "critical", "TF under (2x)"]
for i, z in enumerate(zs):
    u2 = propTF(u1, L1, lam, z)
    ratio = lam * z / (L1 * dx)                      # lam z/(L dx):1 = 临界
    axs[i].imshow(np.cbrt(np.abs(u2) ** 2),
                  extent=[x1[0], x1[-1] + dx, x1[0], x1[-1] + dx],
                  origin="lower", cmap="gray")
    axs[i].set_title(f"z = {z:g} m, {tag[i]}")
    axs[i].set_xlabel("x (m)"); axs[i].set_aspect("equal")
    if i == 1:                                       # z = 1000 m:临界采样档
        u2_ir = propIR(u1, L1, lam, z)
        same = np.max(np.abs(np.abs(u2) ** 2 - np.abs(u2_ir) ** 2)) \
            / np.max(np.abs(u2) ** 2)
        print(f"(d) z=1000(临界) TF 与 IR 辐照度最大相对差 = {same:.1e}")
    print(f"    z={z:g} m: lam z/(L dx) = {ratio:g}")
axs[0].set_ylabel("y (m)")
fig.suptitle(r"Ex 5.1(d): circular aperture Fresnel TF propagation, $I_2^{1/3}$")
fig.tight_layout()
fig.savefig(os.path.join(outdir, "out_ex5_1.png"), dpi=110)
plt.close(fig)

# ======================================================================
# 习题 5.2 三个孔径的临界采样距离(纯计算,式 5.14:z = L1^2/(M lam))
# ======================================================================
print("=" * 60)
print("习题 5.2 临界采样距离 z = L1^2/(M lam),M = 500、lambda = 0.633 um")
M2, lam2 = 500, 0.633e-6
print(f"  (a) L1 = 2 mm   -> z_crit = {(2e-3)**2/(M2*lam2)*100:.2f} cm"
      f"(题给 0.5/1/5 cm 都不是)")
print(f"  (b)(c) L1 = 2.5 cm -> z_crit = {(2.5e-2)**2/(M2*lam2):.3f} m"
      f"(题给 0.5/2/5 m:0.5 在临界之下、2 m 恰近临界、5 m 超过)")

# ======================================================================
# 习题 5.3 功率守恒诊断(数值同 ch5_sampling.py,这里汇总成表)
# ======================================================================
print("=" * 60)
print("习题 5.3 功率诊断(sqr_beam 参数:L = 0.5 m、M = 250、w = 0.051 m)")
Lb, Mb, lam_b = 0.5, 250, 0.5e-6
dxb = Lb / Mb
xb, _ = coords(Lb, Mb)
Xb, Yb = np.meshgrid(xb, xb)
ub = np.where((np.abs(Xb) < 0.051) & (np.abs(Yb) < 0.051), 1.0, 0.0)
P1b = np.sum(np.abs(ub) ** 2) * dxb ** 2
print("     z (m)   P_TF/P1   P_IR/P1")
for z in [1000.0, 2000.0, 4000.0, 20000.0]:
    ptf = np.sum(np.abs(propTF(ub, Lb, lam_b, z)) ** 2) * dxb ** 2 / P1b
    pir = np.sum(np.abs(propIR(ub, Lb, lam_b, z)) ** 2) * dxb ** 2 / P1b
    print(f"    {z:7.0f}  {ptf:.5f}   {pir:.5f}")

# ======================================================================
# 习题 5.4(RS 传播器)提示——笔记给提示,不给数值:
#   (a) H_RS = exp{jkz sqrt(1-(lam fX)^2-(lam fY)^2)}(式 4.19,传播波区),
#       搬进 propTF 框架;h = z/(j lam) * exp(jkr)/r^2(式 4.15)照 propIR 骨架。
#   (b)(c) sqr_beam 衍射角 ~ lam/w ~ 10 urad,全程极傍轴,
#       菲涅耳 vs RS 差别应在数值噪声级;NF = 1.3/2.6/5.2/26。
# ======================================================================

# ======================================================================
# 习题 5.5 高斯光束:w0 = 1 mm、lambda = 0.633 um、L = 15 mm、M = 250
# ======================================================================
print("=" * 60)
print("习题 5.5 高斯光束")
w0 = 1e-3                                            # 束腰 e^-2 半径 (m)
lam = 0.633e-6                                       # 波长 (m)
Lg, Mg = 0.015, 250                                  # 边长 (m)、样本数
dxg = Lg / Mg
xg, _ = coords(Lg, Mg)
Xg, Yg = np.meshgrid(xg, xg)
u0 = np.exp(-(Xg ** 2 + Yg ** 2) / w0 ** 2)          # 式 (5.29),A0 = 1 V/m
zR = np.pi * w0 ** 2 / lam                           # 式 (5.28):瑞利距离
z_crit = Lg ** 2 / (Mg * lam)                        # 式 (5.14):临界采样距离
print(f"zR = pi w0^2/lam = {zR:.3f} m ; 临界采样距离 z = L^2/(M lam) = {z_crit:.2f} m")

fig, axs = plt.subplots(1, 4, figsize=(19.0, 4.6))
row = Mg // 2
for i, z in enumerate([1.0, 5.0, 10.0]):
    u2 = propTF(u0, Lg, lam, z)
    I_num = np.abs(u2[row, :]) ** 2
    wz = w0 * np.sqrt(1 + (z / zR) ** 2)             # 式 (5.27):光束半径
    I_ana = (w0 / wz) ** 2 * np.exp(-2 * xg ** 2 / wz ** 2)   # 式 (5.26)
    rel = np.max(np.abs(I_num - I_ana)) / np.max(I_ana)
    print(f"(a) z = {z:g} m: w(z) = {wz*1e3:.3f} mm, "
          f"数值 vs 解析最大相对偏差 = {rel:.1e}")
    axs[i].plot(xg, I_ana, "--", label="analytic (5.26)")
    axs[i].plot(xg, I_num, "-", label="propTF")
    axs[i].set_title(f"z = {z:g} m, w(z) = {wz*1e3:.2f} mm")
    axs[i].set_xlabel("x (m)"); axs[i].legend(fontsize=8); axs[i].grid(alpha=0.3)
axs[0].set_ylabel("Irradiance")
axs[3].imshow(np.abs(propTF(u0, Lg, lam, 10.0)) ** 2,
              extent=[xg[0], xg[-1] + dxg, xg[0], xg[-1] + dxg],
              origin="lower", cmap="gray")
axs[3].set_title("z = 10 m image"); axs[3].set_xlabel("x (m)")
axs[3].set_aspect("equal")
fig.tight_layout()
fig.savefig(os.path.join(outdir, "out_ex5_5.png"), dpi=110)
plt.close(fig)

# (b) z = 10 m 的源带宽判据:可用带宽 L/(2 lam z);高斯谱 99% 功率带宽
f_avail = Lg / (2 * lam * 10.0)                      # 长距离 regime 的源带宽限制
B99 = erfinv(0.99) / (np.sqrt(2) * np.pi * w0)       # 高斯谱 erf 形式的 99% 带宽
print(f"(b) z = 10 m 可用带宽 L/(2 lam z) = {f_avail:.0f} cyc/m ; "
      f"高斯 99% 功率带宽 = {B99:.0f} cyc/m -> 判据满足")
# (c) 夫琅禾费:I2 = (pi w0^2/(lam z))^2 exp[-2 pi^2 w0^2 r^2/(lam z)^2],
#     与式 (5.26) 的远场渐近 w(z) -> lam z/(pi w0) 一致(推导见笔记)。

# ======================================================================
# 习题 5.6 检验式 (5.15):split-step 10 x 2000 m vs 单次 20000 m
# ======================================================================
print("=" * 60)
print("习题 5.6 split-step 检验式 (5.15)(sqr_beam 参数,每步临界 2000 m)")
Lb2, Mb2 = 0.5, 250
xb2, _ = coords(Lb2, Mb2)
Xb2, Yb2 = np.meshgrid(xb2, xb2)
ub2 = np.where((np.abs(Xb2) < 0.051) & (np.abs(Yb2) < 0.051), 1.0, 0.0)
u_single = propTF(ub2, Lb2, lam_b, 20000.0)
u_split = ub2.astype(complex)
for _ in range(10):                                  # 连续调用 10 次,每步 2000 m
    u_split = propTF(u_split, Lb2, lam_b, 2000.0)
d = np.max(np.abs(u_split) ** 2 - np.abs(u_single) ** 2)
print(f"单次 vs 拆步辐照度最大差 = {d:.1e}(浮点误差级 -> 结果相同)")

# ======================================================================
# 习题 5.8 方孔 98% 谱功率判据(纯逻辑,串三根线,答案:约 200/256)
# ======================================================================
print("=" * 60)
print("习题 5.8:98% 谱功率 -> B1 = 5.98/w;夫琅禾费相位采样 -> dx1 >= lam z/L1;")
print("  夫琅禾费 regime -> lam z >> w^2。判据 B1 < 1/(2 dx1) -> w/dx1 > 11.96,")
print("  代入 M >> (w/dx1)^2(式 5.24)得 M >> 143,取整 M ~ 200(2 的幂取 256)。")

print("=" * 60)
print("已保存 out_ex5_1.png / out_ex5_5.png")
