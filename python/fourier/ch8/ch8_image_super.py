# -*- coding: utf-8 -*-
"""书 8.7.2 节 image_super:空间变非相干成像的叠加仿真。

理想像 Ig 的每个样本 (u0,v0) 处:生成该像点的 PSF(式 8.7) -> 归一化体积
为 1(等价 OTF/MTF 归一化,书式 7.25) -> circshift 平移到像点位置 -> 乘
Ig(u0,v0) 后累加(式 8.17 的离散实现)。像差参数与 psf_map 相同
(0.5λ 球差 + 1λ 彗差 + 1.5λ 像散,f/20 出瞳)。

原书读入随书代码分发的 USAF1951B250.png;Python 侧用合成 USAF 风格靶
(make_usaf_target)替代。全书脚本默认 M=250,62,500 个 PSF,书作者笔记本
约 1 小时;命令行可传 M(如 `python ch8_image_super.py 100`)做快速验证。

MATLAB/Python 差异:
  - 书第 36-39 行用循环索引 n、m 推出 v0、u0(MATLAB 行列序:行 = v);
    Python 用 nj、mi(0-based)同构:u0 = (mi - M//2)/(M//2)。
  - np.roll 轴 (0,1) = (v 行, u 列),与 MATLAB circshift [行,列] 一致。
运行产物:out_fig8_imagesuper.png + 控制台验证。
"""
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, HERE)
from base import coords, circle                          # noqa: E402
from ch8_aberrations import seidel_5, psf_from_H, make_usaf_target  # noqa: E402

M = int(sys.argv[1]) if len(sys.argv) > 1 else 250       # 书默认 250
L = 1e-3                # image plane side length
du = L / M              # sample interval
u, _ = coords(L, M)     # coordinates
v = u

lam = 0.5e-6            # wavelength
k = 2 * np.pi / lam     # wavenumber
wxp = 2.5e-3            # exit pupil radius
zxp = 100e-3            # exit pupil distance
fnum = zxp / (2 * wxp)  # exit pupil f-number
twof0 = 1 / (lam * fnum)    # inc cutoff freq
fN = 1 / (2 * du)           # Nyquist frequency
print(f"M = {M}, fnum = {fnum:.0f}, 2f0 = {twof0:.4g}, fN = {fN:.4g}, "
      f"式(8.6): {twof0 <= fN}")

coeffs = (0 * lam, 0.5 * lam, 1 * lam, 1.5 * lam, 0 * lam, 0 * lam)

# 理想像:合成 USAF 风格靶(原书 USAF1951B250.png,行 0 = 图顶 ->
# flipud;合成靶已按"行 0 = 底"画出,直接可用)
Ig = make_usaf_target(M)
Ig = Ig / Ig.max()                                       # 书 5 行:normalize

# 频率坐标预移位(书 30-31 行,省去循环内 fftshift)
_, fu = coords(L, M)
fu = np.fft.fftshift(fu)
FU, FV = np.meshgrid(fu, fu)

I = np.zeros((M, M))                                     # 书 34 行
t0 = time.time()
h2 = np.empty((M, M))                                    # 复用缓冲,少一次分配
for nj in range(M):                                      # 书 36 行:行 = v 方向
    v0 = (nj - M // 2) / (M // 2)                        # 书 37 行
    for mi in range(M):                                  # 书 38 行:列 = u 方向
        u0 = (mi - M // 2) / (M // 2)                    # 书 39 行
        # wavefront + coherent TF + PSF(书 41-48 行)
        W = seidel_5(u0, v0, -2 * lam * fnum * FU, -2 * lam * fnum * FV, *coeffs)
        H = circle(FU, FV, wxp / (lam * zxp)) * np.exp(-1j * k * W)
        np.abs(np.fft.ifft2(H), out=h2)                  # h = |IFFT{H}|(FFT 序)
        h2 *= h2                                         # 书 48 行:h2 = |h|^2
        h2 /= h2.sum()                                   # 书 49 行:体积归一到 1
        h2c = np.fft.ifftshift(h2)                       # 回中(自然序 -> 居中)
        # shift 到像点位置(书 51 行;roll 量 = 索引差)
        h2c = np.roll(h2c, (nj - M // 2, mi - M // 2), axis=(0, 1))
        I += Ig[nj, mi] * h2c                            # 书 53 行:叠加
dt = time.time() - t0
print(f"叠加 {M*M} 个 PSF 耗时 {dt:.1f} s"
      f"(书默认 M=250 约 62500 个 PSF,书作者笔记本 ~1 小时)")

fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))
ext = [u[0] * 1e3, u[-1] * 1e3, v[0] * 1e3, v[-1] * 1e3]   # mm
axes[0].imshow(Ig, origin="lower", extent=ext, cmap="gray")
axes[0].set_title("Ideal image (synthetic USAF-style target)")
axes[1].imshow(I ** (1 / 3), origin="lower", extent=ext, cmap="gray")
axes[1].set_title(f"Superposition image (M={M}, nthroot 3)")
for ax in axes:
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig8_imagesuper.png"), dpi=150)
print("saved out_fig8_imagesuper.png")

# 验证:总能量守恒(PSF 体积归一 -> 叠加像总功率 = 理想像总功率)
print(f"理想像总功率 {Ig.sum():.6f},叠加像总功率 {I.sum():.6f},"
      f"比值 {I.sum()/Ig.sum():.6f}(应=1)")
