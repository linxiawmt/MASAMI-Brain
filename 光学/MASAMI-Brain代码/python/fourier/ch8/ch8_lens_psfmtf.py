# -*- coding: utf-8 -*-
"""书 8.5 节 lens_psfmtt:f/5 平凸透镜(Newport KPX094)的 PSF 与 MTF。

原书脚本逐行对应见 matlab/fourier/ch8/lens_psfmtf.m;本脚本把"轴向像点
(0,0) / 最大像高点 (0,1)"两种情形各跑一遍,输出:
  - 相干传递函数相位(式 8.3/8.4,查式 8.6 采样判据:瞳直径 <= 阵列之半)
  - PSF(式 8.7)图样与 u/v 剖面
  - MTF(式 8.9)与一维衍射限解析曲线(式 7.31)对比
像差系数来自 ZEMAX(书表 8.2,单位:波长数 x 0.55 um)。

MATLAB/Python 差异:
  - MATLAB 剖面取 h2(M/2+1,:)与 h2(:,M/2+1);numpy 取 h2[M//2,:]与 h2[:,M//2]。
  - MATLAB imagesc(u,v,...)+axis xy 对应 Python imshow(...,origin='lower',
    extent=[u0,u1,v0,v1])(数组行 = v)。
运行产物:out_fig8_lens.png(本章目录)+ 控制台数值验证。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, HERE)            # python -I 隔离模式不含脚本目录,显式加入
from base import coords, circle, rect                      # noqa: E402
from ch8_aberrations import (seidel_5, psf_from_H,         # noqa: E402
                             mtf_from_psf, diffraction_limited_mtf_1d)

HERE = os.path.dirname(os.path.abspath(__file__))

# ----- 参数(书 lens_psfmtf 第 6-27 行) -----
M = 1024                # sample #
L = 1e-3                # image plane side length
du = L / M              # sample interval
u, _ = coords(L, M)     # coordinates(-L/2 .. L/2-du,居中)
v = u

lam = 0.55e-6           # wavelength
k = 2 * np.pi / lam     # wavenumber
Dxp = 20e-3
wxp = Dxp / 2           # exit pupil radius
zxp = 100e-3            # exit pupil distance
fnum = zxp / (2 * wxp)  # exit pupil f-number (= 5)
lz = lam * zxp
twof0 = 1 / (lam * fnum)    # incoh cutoff freq
fN = 1 / (2 * du)           # Nyquist frequency

# ZEMAX Seidel 系数(书表 8.2;"波长数" x lambda)
wd = 0 * lam
w040 = 4.963 * lam
w131 = 2.637 * lam
w222 = 9.025 * lam
w220 = 7.536 * lam
w311 = 0.157 * lam
coeffs = (wd, w040, w131, w222, w220, w311)

# 频率坐标(居中排列;lens_psfmtf 第 29-30 行)
# 注意:base.coords 返回 (空间坐标 x, 频率坐标 f),频率轴是第二个输出!
_, fu = coords(L, M)        # -1/(2du) .. 1/(2du)-1/L,步距 1/L
FU, FV = np.meshgrid(fu, fu)


def run_point(u0, v0):
    """对归一化像点 (u0,v0) 返回 H, h2, MTF(书第 33-56 行)。"""
    # wavefront:归一化瞳坐标 xhat = -lz*Fu/wxp = -2*lambda*fnum*Fu
    W = seidel_5(u0, v0, -lz * FU / wxp, -lz * FV / wxp, *coeffs)   # 书 33-34 行
    # coherent transfer function(式 8.3/8.4):circ(f_r/f0)·exp(-jkW)
    H = circle(FU, FV, wxp / lz) * np.exp(-1j * k * W)              # 书 37 行
    # PSF(式 8.7):H 为居中排列 -> ifftshift -> ifft2 -> fftshift 回中
    h2 = psf_from_H(H, pre_shifted=False)                            # 书 43 行
    MTF = mtf_from_psf(h2)                                           # 书 54-56 行
    return H, h2, MTF


def practical_cutoff(freq_axis, mtf_1d, level=0.05):
    """MTF 实用截止频率:从 0 频起找 MTF 首次跌破 level 的频率。"""
    half = freq_axis >= 0
    f = freq_axis[half]
    m = mtf_1d[half]
    below = np.where(m < level)[0]
    return f[below[0]] / 1e3 if below.size else np.nan   # cyc/mm


results = {}
for tag, (u0, v0) in {"axial": (0, 0), "field": (0, 1)}.items():
    H, h2, MTF = run_point(u0, v0)
    results[tag] = (H, h2, MTF)
    peak = h2.max()
    cut_u = practical_cutoff(fu, MTF[M // 2, :])
    cut_v = practical_cutoff(fu, MTF[:, M // 2])
    cut_u2 = practical_cutoff(fu, MTF[M // 2, :], level=0.02)
    cut_v2 = practical_cutoff(fu, MTF[:, M // 2], level=0.02)
    print(f"[{tag}] (u0,v0)=({u0},{v0})  PSF peak = {peak:.4g}")
    print(f"    MTF 实用截止(0.05 电平): u 轴 {cut_u:.1f} cyc/mm, v 轴 {cut_v:.1f} cyc/mm")
    print(f"    MTF 实用截止(0.02 电平): u 轴 {cut_u2:.1f} cyc/mm, v 轴 {cut_v2:.1f} cyc/mm"
          f"(书正文目测:轴上 ~120-130,最大像高 u ~20 / v ~10 cyc/mm)")
    for fq in (10, 20, 120, 130):
        iu = np.argmin(np.abs(fu - fq * 1e3))
        print(f"    MTF @ {fq} cyc/mm: u 轴 {MTF[M//2, iu]:.3f}, v 轴 {MTF[iu, M//2]:.3f}")

# ----- 数值验证输出 -----
print(f"2f0 = {twof0:.4g} cyc/m = {twof0/1e3:.1f} cyc/mm;  fN = {fN:.4g} cyc/m")
print(f"式(8.6) 判据 2f0 <= fN : {twof0 <= fN};瞳直径(2f0)占阵列频率宽度比例 = "
      f"{twof0/(1/du):.3f}(< 0.5 满足)")
# 习题 8.8(a):无像差 PSF 峰值(轴向)
H0 = circle(FU, FV, wxp / lz)
h2_unab = psf_from_H(H0, pre_shifted=False)
S_axial = results["axial"][1].max() / h2_unab.max()
print(f"习题 8.8(e):无像差峰值 Iunab = {h2_unab.max():.4g},"
      f"轴上 Strehl S = {S_axial:.3g}")
# 衍射限 MTF 与解析式逐点对比(验证 MTF 口径)
h2_dl = psf_from_H(H0, pre_shifted=False)
MTF_dl = mtf_from_psf(h2_dl)
MTF_an = diffraction_limited_mtf_1d(fu, twof0)
rel = np.max(np.abs(MTF_dl[M // 2, fu >= 0] - MTF_an[fu >= 0]))
print(f"衍射限 MTF(u 轴)vs 式(7.31) 解析:最大偏差 = {rel:.2e}")

# ----- 出图(2 x 3:左列轴向,右列最大像高;图内标签用英文,避免字体缺字) -----
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
ext = [u[0] * 1e3, u[-1] * 1e3, v[0] * 1e3, v[-1] * 1e3]   # mm
for col, tag in enumerate(["axial", "field"]):
    H, h2, MTF = results[tag]
    im = axes[0, col].imshow(np.sqrt(h2), origin="lower", extent=ext,
                             cmap="gray", vmin=0, vmax=np.sqrt(h2).max() * 0.05)
    axes[0, col].set_title("PSF (sqrt stretch, "
                           + ("on-axis" if col == 0 else "max field") + ")")
    axes[0, col].set_xlabel("u (mm)"); axes[0, col].set_ylabel("v (mm)")
    plt.colorbar(im, ax=axes[0, col], fraction=0.046)
    fhalf = fu >= 0
    axes[1, col].plot(fu[fhalf] / 1e3, MTF[M // 2, fhalf], label="u MTF")
    axes[1, col].plot(fu[fhalf] / 1e3, MTF[:, M // 2][fhalf], "--", label="v MTF")
    axes[1, col].plot(fu[fhalf] / 1e3, MTF_an[fhalf], ":", label="diff. limit")
    axes[1, col].set_xlim(0, 150); axes[1, col].set_ylim(0, 1)
    axes[1, col].set_xlabel("f (cyc/mm)"); axes[1, col].set_ylabel("Modulation")
    axes[1, col].set_title("MTF"); axes[1, col].legend(fontsize=8)
# 相位显示:瞳外 |H|=0,但 0*exp(jθ) 会留下"带符号零"(如 -0.0-0.0j),
# np.angle/angle 对其返回 ±π,造成书 Fig.8.5(a) 所说的"display artifact"。
# 显示前把 |H|=0 处的相位置 0(书 8.5 节评注 c:瞳外冗余值"can be removed later")。
angH = np.angle(results["axial"][0])
angH[results["axial"][0] == 0] = 0.0
axes[0, 2].imshow(angH, origin="lower", extent=ext, cmap="gray")
axes[0, 2].set_title("Coherent TF phase (on-axis)")
axes[1, 2].plot(u * 1e3, results["field"][1][:, M // 2], label="v profile")
axes[1, 2].plot(u * 1e3, results["field"][1][M // 2, :], label="u profile")
axes[1, 2].set_title("Max-field PSF profiles"); axes[1, 2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig8_lens.png"), dpi=150)
print("saved out_fig8_lens.png")
