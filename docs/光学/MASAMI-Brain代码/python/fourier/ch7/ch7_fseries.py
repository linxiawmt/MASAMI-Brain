# -*- coding: utf-8 -*-
"""书图 7.9(相干)与图 7.14(非相干):f/5、f/10、f/20 系列的谱与像。

理想图像的奈奎斯特频率 fN = 4.17e5 cyc/m 固定不变(L=0.3 mm、M=250),
改变出瞳半径 wxp(即 f/#)只改变传递函数的截止:
  相干:  f0 = 1/(2*lam*(f/#))            —— f/5: 2e5, f/10: 1e5, f/20: 5e4
  非相干: 2f0 = 1/(lam*(f/#))            —— f/5: 4e5, f/10: 2e5, f/20: 1e5
谱幅度用对数(相干)/五次方根(非相干,书图 7.14)显示以突出低对比细节。
输出:out_fig7_9.png、out_fig7_14.png。
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
from ch7_imaging import (freq_coords, H_circ, otf_from_H,  # noqa: E402
                         coh_image, incoh_image)
from ch7_usaf_chart import ensure_chart, L_IMG, N_PX     # noqa: E402

Ig = ensure_chart()
Ig = np.flipud(Ig)
Ig = Ig / Ig.max()

L, M = L_IMG, N_PX
du = L / M
u, _ = coords(L, M)
lam = 0.5e-6
fu = freq_coords(L, M)
fN = 1.0 / (2 * du)

cases = [(5.0, 2.0e5), (10.0, 1.0e5), (20.0, 0.5e5)]   # (f/#, 相干截止 f0)

# ---------- 相干系列(书图 7.9:谱 log 幅度 | 像) ----------
fig, axs = plt.subplots(4, 2, figsize=(8.6, 12.0))
FI_ideal = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Ig)))
log_ideal = np.log10(np.abs(FI_ideal) / np.abs(FI_ideal).max() + 1e-12)
axs[0, 0].imshow(log_ideal, vmin=-6, vmax=0,
                 extent=[fu[0], fu[-1], fu[0], fu[-1]], origin="lower", cmap="gray")
axs[0, 0].set_title(f"Ideal spectrum  fN = {fN:.2e}")
axs[0, 1].imshow(Ig, extent=[u[0], u[-1] + du, u[0], u[-1] + du],
                 origin="lower", cmap="gray")
axs[0, 1].set_title("Ideal image")
for k, (fn, f0) in enumerate(cases, start=1):
    H = H_circ(fu, f0)
    ui, Gi = coh_image(Ig ** 0.5, H)         # 理想场 = sqrt(Ig)(二值图同 Ig)
    Ii = np.abs(ui) ** 2
    log_sp = np.log10(np.abs(Gi) / np.abs(Gi).max() + 1e-12)
    axs[k, 0].imshow(np.where(H > 0, log_sp, -8), vmin=-6, vmax=0,
                     extent=[fu[0], fu[-1], fu[0], fu[-1]],
                     origin="lower", cmap="gray")
    axs[k, 0].set_title(f"f/{fn:g}   f0 = {f0:.1e}")
    axs[k, 1].imshow(np.sqrt(Ii), extent=[u[0], u[-1] + du, u[0], u[-1] + du],
                     origin="lower", cmap="gray", vmin=0, vmax=1)
    axs[k, 1].set_title(f"f/{fn:g} coherent image")
for ax in axs.ravel():
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_9.png"), dpi=100)
plt.close(fig)

# ---------- 非相干系列(书图 7.14:谱五次方根 | 像) ----------
fig, axs = plt.subplots(4, 2, figsize=(8.6, 12.0))
FI = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Ig)))
axs[0, 0].imshow((np.abs(FI) / np.abs(FI).max()) ** 0.2,
                 extent=[fu[0], fu[-1], fu[0], fu[-1]], origin="lower", cmap="gray")
axs[0, 0].set_title(f"Ideal spectrum  fN = {fN:.2e}")
axs[0, 1].imshow(Ig, extent=[u[0], u[-1] + du, u[0], u[-1] + du],
                 origin="lower", cmap="gray")
axs[0, 1].set_title("Ideal image")
for k, (fn, f0) in enumerate(cases, start=1):
    OTF = otf_from_H(H_circ(fu, f0))
    Ii = incoh_image(Ig, OTF)
    Gi = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Ii)))
    axs[k, 0].imshow((np.abs(Gi) / np.abs(Gi).max()) ** 0.2,
                     extent=[fu[0], fu[-1], fu[0], fu[-1]],
                     origin="lower", cmap="gray")
    axs[k, 0].set_title(f"f/{fn:g}   2f0 = {2*f0:.1e}")
    axs[k, 1].imshow(np.sqrt(Ii), extent=[u[0], u[-1] + du, u[0], u[-1] + du],
                     origin="lower", cmap="gray", vmin=0, vmax=1)
    axs[k, 1].set_title(f"f/{fn:g} incoherent image")
for ax in axs.ravel():
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_14.png"), dpi=100)
plt.close(fig)

print("== 书图 7.9 / 7.14 系列 ==")
print(f"fN = {fN:.3e} cyc/m(理想图像不变)")
for fn, f0 in cases:
    print(f"f/{fn:g}: 相干截止 f0 = {f0:.2e}, 非相干截止 2f0 = {2*f0:.2e} cyc/m")
print("已保存 out_fig7_9.png / out_fig7_14.png")
