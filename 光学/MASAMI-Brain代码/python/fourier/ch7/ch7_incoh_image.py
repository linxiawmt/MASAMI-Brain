# -*- coding: utf-8 -*-
"""书 7.3.3 节 incoh_image 非相干成像例子的 Python 复现(书图 7.11/7.12/7.13)
+ 式 (7.29)/(7.31) OTF 解析式对拍 + 式 (7.17) "免缠绕"条件的数值演示。

参数与 7.2.3 节 coh_image 完全相同(f/10、L=0.3 mm、M=250)。
输出:out_fig7_11.png(OTF 曲面)、out_fig7_12.png(模拟非相干像+剖面)、
     out_fig7_otfcheck.png(OTF 解析 vs FFT;f/10 合格 vs f/4 违反判据 7.17
     时的环绕污染)+ 命令行数值检查。

MATLAB/Python 差异:
  - 书第 22-23 行的 OTF 留在 FFT 排列,与 fft2(fftshift(Ig)) 对乘;
    本模块 otf_from_H/incoh_image 原样保留该约定(见 ch7_imaging.py)。
  - 书第 36 行 Ii=real(Ii); mask=Ii>=0; Ii=mask.*Ii —— numpy 侧在
    incoh_image() 里等价实现。
  - MATLAB 的 OTF(1,1)(1-based)= numpy OTF[0,0](0-based),都是 FFT
    排列的零频样本。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
sys.path.insert(0, os.path.join(_here, "..", "common"))
from base import coords                                  # noqa: E402
from ch7_imaging import (freq_coords, H_circ, H_square, otf_from_H,  # noqa: E402
                         incoh_image, otf_circ_analytic, otf_square_analytic)
from ch7_usaf_chart import ensure_chart, L_IMG, N_PX     # noqa: E402

# ---------- 书第 3-16 行:理想图像与系统参数(与 coh_image 相同) ----------
Ig = ensure_chart()
Ig = np.flipud(Ig)
Ig = Ig / Ig.max()

L = L_IMG
M = N_PX
du = L / M
u, _ = coords(L, M)
v = u

lam = 0.5e-6
wxp = 6.25e-3
zxp = 125e-3
f0 = wxp / (lam * zxp)               # 相干截止 1e5 cyc/m(书 16 行)
fu = freq_coords(L, M)
fN = 1.0 / (2 * du)

# ---------- 书第 18-23 行:相干传递函数 -> OTF ----------
H = H_circ(fu, f0)                   # 书 21 行
OTF = otf_from_H(H)                  # 书 22-23 行(FFT 排列、零频归一)

# ---------- 书第 32-36 行:式 (7.33) 非相干成像 ----------
Ii = incoh_image(Ig, OTF)

row = M // 2 + int(round(+0.2e-4 / du))   # 书 45-46 行:v=+0.2e-4 的行号

# ================= 出图(书图 7.11/7.12/7.13) =================
fig = plt.figure(figsize=(5.6, 4.6))
ax = fig.add_subplot(projection="3d")
FU, FV = np.meshgrid(fu, fu)
ax.plot_surface(FU, FV, np.fft.fftshift(np.abs(OTF)), cmap="gray",
                rstride=2, cstride=2, linewidth=0, antialiased=True)
ax.set_xlabel("fv (cyc/m)"); ax.set_ylabel("fu (cyc/m)")
ax.set_title("OTF (diffraction-limited, circular pupil)")
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_11.png"), dpi=110)
plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.5))
axs[0].imshow(np.sqrt(Ii), extent=[u[0], u[-1] + du, v[0], v[-1] + du],
              origin="lower", cmap="gray", vmin=0, vmax=1)
axs[0].set_xlabel("u (m)"); axs[0].set_ylabel("v (m)")
axs[0].set_title("simulated incoherent image (sqrt Ii)")
axs[0].annotate("Group -1\nElement 3", xy=(-0.30e-4, 0.2e-4),
                xytext=(0.0e-4, 0.85e-4), color="w", fontsize=8,
                arrowprops=dict(arrowstyle="->", color="w"))
axs[1].plot(u, Ii[row, :], label="simulated")
axs[1].plot(u, Ig[row, :], ":", label="ideal")
axs[1].set_xlabel("u (m)"); axs[1].set_ylabel("Irradiance")
axs[1].set_title("row profile at v = +0.2e-4 m")
axs[1].legend(); axs[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_12.png"), dpi=110)
plt.close(fig)

# ================= OTF 解析对拍(式 7.31 圆瞳 / 式 7.29 方瞳) =================
OTF_c = np.fft.fftshift(OTF)          # 转回居中排列再与解析式比
otf_a = otf_circ_analytic(fu, f0)
err_circ = np.abs(OTF_c - otf_a).max()

Hsq = H_square(fu, f0)                # 方瞳(式 7.9)
OTFsq_c = np.fft.fftshift(otf_from_H(Hsq))
otf_sq = otf_square_analytic(fu, f0)
err_sq = np.abs(OTFsq_c - otf_sq).max()

# OTF 支撑半径(应 = 非相干截止 2f0)
fr = np.hypot(*np.meshgrid(fu, fu))
sup = fr[OTF_c > 1e-6].max()

# ================= 式 (7.17) 的第二张面孔:FFT 自相关免缠绕 =================
# OTF 用 FFT 自相关计算 = 周期延拓下的循环自相关;相邻副本间距 1/du,
# 免缠绕条件 2f0 <= 1/du - 2f0,即 2f0 <= fN —— 恰为判据 (7.17)。
# 用 f/4(2f0 = 5e5 > fN = 4.17e5)演示违反时的环绕污染。
H4 = H_circ(fu, 1.0 / (lam * 4.0) / 2)        # f/#=4: 2f0 = 1/(lam*4) = 5e5
OTF4_c = np.fft.fftshift(otf_from_H(H4))
otf4_a = otf_circ_analytic(fu, 1.0 / (lam * 4.0) / 2)
err4 = np.abs(OTF4_c - otf4_a).max()
sl = np.abs(fu) < 3.0e5                       # 沿 fV=0 切片对比

fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2))
axs[0].plot(fu, OTF_c[M // 2, :].real, label="FFT autocorrelation (f/10)")
axs[0].plot(fu, otf_a[M // 2, :], "--", label="analytic Eq.(7.31)")
axs[0].set_title(f"f/10: 2f0=2e5 <= fN=4.17e5, max|diff| = {err_circ:.1e}")
axs[0].set_xlabel("fu (cyc/m)"); axs[0].legend(); axs[0].grid(alpha=0.3)
axs[1].plot(fu[sl], OTF4_c[M // 2, sl].real, label="FFT autocorrelation (f/4)")
axs[1].plot(fu[sl], otf4_a[M // 2, sl], "--", label="analytic Eq.(7.31)")
axs[1].set_title(f"f/4: 2f0=5e5 > fN, wrap-around, max|diff| = {err4:.1e}")
axs[1].set_xlabel("fu (cyc/m)"); axs[1].legend(); axs[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_otfcheck.png"), dpi=110)
plt.close(fig)

# ================= 剖面里组-1/单元3 的可分辨性(书 7.3.3 节判读) =================
# 本板的组-1/单元3 频率 = 2e5*2^(-1+3/6) = 1.414e5 cyc/m(书板判读约 1.7e5)
f_e13 = 2.0e5 * 2.0 ** (-1 + 3.0 / 6.0)

print("== 书 7.3.3 节 incoh_image 数值检查 ==")
print(f"OTF 零频值 OTF(0,0)       = {OTF_c[M//2, M//2]:.6f}  (归一化要求 = 1)")
print(f"OTF 最小值(衍射限应>=0)  = {OTF_c.min():.3e}")
print(f"OTF 支撑半径              = {sup:.3e} cyc/m  (理论 2f0 = {2*f0:.3e})")
print(f"圆瞳 OTF vs 式(7.31)      最大差 = {err_circ:.2e}")
print(f"方瞳 OTF vs 式(7.29)      最大差 = {err_sq:.2e}")
print(f"f/4 违反判据 7.17 时       最大差 = {err4:.2e}(环绕污染,见图右)")
print(f"组-1/单元3 频率           = {f_e13:.3e} cyc/m  (< 2f0 = 2e5,非相干可分辨;"
      f" > f0 = 1e5,相干不可分辨)")
seg = Ii[row, :]
loc = (u > -3.9e-5) & (u < -1.9e-5)          # 组-1 列的单元3 三条竖条范围
seg_e13 = seg[loc]
print(f"组-1/单元3 剖面调制度     = {(seg_e13.max()-seg_e13.min())/seg_e13.max():.2f}"
      f"(理想 1.0;OTF 在 1.41e5 处传函 0.18,条纹仍可辨认)")
loc23 = (u > -1.05e-4) & (u < -6.0e-5)       # 同一行穿过的组-2/单元3(7.07e4)
seg_e23 = seg[loc23]
print(f"组-2/单元3 剖面调制度     = {(seg_e23.max()-seg_e23.min())/seg_e23.max():.2f}"
      f"(OTF 在 7.07e4 处传函 0.56,充分分辨)")
row_e1 = M // 2 + int(round(+1.0e-4 / du))   # 组-2/单元1 大条所在行
seg1 = Ii[row_e1, :][(u > -1.14e-4) & (u < -0.62e-4)]
print(f"组-2/单元1 处最大值       = {seg1.max():.3f}(相干像同处 1.50 —— 非相干无振铃)")
print("已保存 out_fig7_11/12.png / out_fig7_otfcheck.png")
