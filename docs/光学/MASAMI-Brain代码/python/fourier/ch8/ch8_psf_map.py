# -*- coding: utf-8 -*-
"""书 8.7.1 节 psf_map:像面上不同位置的 PSF 分布图(f/20 出瞳)。

对 7x7 个归一化像点 (u0,v0) 各生成一幅 PSF(式 8.7),用 circshift 平移到
像面对应位置后叠加显示。像差夸大取 0.5λ 球差 + 1λ 彗差 + 1.5λ 像散
(书 8.7.1 节原话),以凸显 PSF 随像点位置的变化。

省时技巧(书 8.7.1 评注 b):先把频率坐标 fu=fftshift(fu) 预移位,循环内
H 天然处于 FFT 序,ifft2 前不再需要 fftshift。

MATLAB/Python 差异:
  - MATLAB circshift(h2,[r,c]) 第一个变量是行(= v 方向);numpy
    np.roll(h2, (r, c), axis=(0,1)) 轴含义相同。
  - round 半整数:本题 v0*M/2 会出现 ±87.5,numpy/MATLAB 都给 ±88(偶数),
    结果一致;其他取值建议显式 int(np.floor(x+0.5)) 以避开银行家舍入差异。
运行产物:out_fig8_psfmap.png + 控制台验证。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, HERE)
from base import coords, circle                           # noqa: E402
from ch8_aberrations import seidel_5, psf_from_H          # noqa: E402

# ----- 参数(书 psf_map 第 3-27 行) -----
M = 250
L = 1e-3                # image plane side length
du = L / M              # sample interval
u, _ = coords(L, M)     # coordinates
v = u

lam = 0.5e-6            # wavelength
k = 2 * np.pi / lam     # wavenumber
wxp = 2.5e-3            # exit pupil radius
zxp = 100e-3            # exit pupil distance
fnum = zxp / (2 * wxp)  # = 20
twof0 = 1 / (lam * fnum)    # incoh cutoff freq
fN = 1 / (2 * du)           # Nyquist frequency
print(f"fnum = {fnum:.0f}, 2f0 = {twof0:.4g} 1/m, fN = {fN:.4g} 1/m, "
      f"式(8.6) 2f0<=fN: {twof0 <= fN}")

# 像差系数(书第 18-23 行:0.5λ 球差、1λ 彗差、1.5λ 像散)
coeffs = (0 * lam, 0.5 * lam, 1 * lam, 1.5 * lam, 0 * lam, 0 * lam)

# 频率坐标:预移位(书第 25-26 行 fu=fftshift(fu))
_, fu = coords(L, M)
fu = np.fft.fftshift(fu)
FU, FV = np.meshgrid(fu, fu)

I = np.zeros((M, M))                                    # 书 29 行
for u0 in np.arange(-0.7, 0.7001, 0.7 / 3):             # 书 31 行:7 个位置
    for v0 in np.arange(-0.7, 0.7001, 0.7 / 3):         # 书 32 行
        # wavefront:归一化瞳坐标 xhat = -2*lambda*fnum*Fu(= -lz*Fu/wxp)
        W = seidel_5(u0, v0, -2 * lam * fnum * FU,
                     -2 * lam * fnum * FV, *coeffs)     # 书 34-36 行
        # coherent transfer function(式 8.3/8.4):书 circ(f_r*2*lambda*fnum)
        # = circ(f_r/f0),f0 = wxp/(lam*zxp) = 1/(2*lambda*fnum) —— base.circle
        # 的第三参是"半径",故传 f0(切勿传 1/twof0,那是长度量纲的 λf/#)
        H = circle(FU, FV, wxp / (lam * zxp)) * np.exp(-1j * k * W)   # 书 38-39 行
        # PSF:H 已是 FFT 序 -> 直接 ifft2,ifftshift 回中
        h2 = psf_from_H(H, pre_shifted=True)            # 书 41 行
        # shift PSF to image plane position(书 43-44 行)
        h2 = np.roll(h2, (int(np.round(v0 * M / 2)),
                          int(np.round(u0 * M / 2))), axis=(0, 1))
        I = I + h2                                      # 书 46 行:仅摆放,非卷积!

# ----- 出图(书图 8.9) -----
fig, ax = plt.subplots(figsize=(6.4, 5.4))
ext = [u[0] * 1e3, u[-1] * 1e3, v[0] * 1e3, v[-1] * 1e3]   # mm
im = ax.imshow(I ** 0.25, origin="lower", extent=ext, cmap="gray")
ax.set_xlabel("u (mm)"); ax.set_ylabel("v (mm)")
ax.set_title("PSF map (nthroot 4): 0.5λ spherical, 1λ coma, 1.5λ astig")
fig.colorbar(im, ax=ax, fraction=0.046)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig8_psfmap.png"), dpi=150)
print("saved out_fig8_psfmap.png")

# 验证:中心 PSF 对称(0 频像点只受球差影响,书 8.2 节评注 a),
# 边缘 PSF 沿 v 拉长(彗差 + 像散随 |u0| 增大)
c = I.shape[0] // 2
cen = I[c - 6:c + 7, c - 6:c + 7]
print(f"中心 PSF 峰值 {I[c, c]:.4g};中心 13x13 邻域左右对称差 = "
      f"{np.max(np.abs(cen - cen[:, ::-1])) / I[c, c]:.2e}(球差为主应近对称)")
