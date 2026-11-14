# -*- coding: utf-8 -*-
"""书 6.5 节习题 6.2/6.3/6.5/6.7/6.8:Python 实跑验证。

6.2  tilt + focus 同时作用(sqr_beam,TF,z=2000 m)
6.3  柱面聚焦(焦线)+ (c) 焦面不是二维夫琅禾费的说明
6.5  波带片(式 6.29 cos 型 / 式 6.30 二值型),z=5/8/10/12 m,
     与同参数理想正透镜(f=10 m、半径 w)的焦面对比
6.7  相位光栅(式 6.7(a)),m=pi、pi/2、pi/4 的零阶功率 J0(m)^2
6.8  三角波光栅(式 6.32)一维模型 vs 解析(阶权重 ∝ sinc²(n/2))
输出:out_fig6_ex63.png、out_fig6_ex65.png、out_fig6_ex68.png
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
from base import coords, sinc, circle          # noqa: E402
from ch5_prop import propTF, propFF            # noqa: E402
from ch6_transmittance import (tilt, focus, focus_cyl, zone_plate,  # noqa: E402
                               tri_wave_1d)
from scipy.special import j0, jv

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = np.pi / 180


def lens_ta_plane(X, Y, lam, f, wL):
    """习题 6.5(h) 的理想正透镜透过率(同 f、同口径,式 6.12)。"""
    k = 2 * np.pi / lam
    return np.exp(-1j * k / (2 * f) * (X ** 2 + Y ** 2)) * circle(X, Y, wL)


def fwhm(prof, x):
    """半高全宽(以峰为基准的简单估计)。"""
    pk = prof.max()
    above = x[prof > pk / 2]
    return above[-1] - above[0]

# ================= 习题 6.2:tilt + focus =================
L1 = 0.5; M = 250; dx1 = L1 / M
x1, _ = coords(L1, M)
lam = 0.5e-6; w = 0.051; z = 2000.0
X1, Y1 = np.meshgrid(x1, x1)
u1 = np.where((np.abs(X1) < w) & (np.abs(Y1) < w), 1.0, 0.0)

u1t = tilt(u1, L1, lam, 5.0e-5, 45 * DEG)
u1tf = focus(u1t, L1, lam, 4000.0)             # zf=4000 m(习题指定)
u2 = propTF(u1tf, L1, lam, z)
I2 = np.abs(u2) ** 2
P2 = I2.sum()
xc = (I2.sum(axis=0) * x1).sum() / P2
yc = (I2.sum(axis=1) * x1).sum() / P2
I2_pure = np.abs(propTF(focus(u1, L1, lam, 4000.0), L1, lam, z)) ** 2
print(f"[6.2] 质心 = ({xc*1e3:.2f}, {yc*1e3:.2f}) mm(理论 45° 方向 "
      f"{z*np.tan(5e-5)*1e3:.1f} mm → 分量 {z*np.tan(5e-5)*np.cos(45*DEG)*1e3:.2f});"
      f"峰值 tilt+focus/zf=4000 = {I2.max():.1f} vs 纯 focus = {I2_pure.max():.1f}")

# ================= 习题 6.3:柱面聚焦 =================
u1c = focus_cyl(u1, L1, lam, 2000.0)
I2c = np.abs(propTF(u1c, L1, lam, 2000.0)) ** 2
profc = M // 2
fx_width = x1[I2c[profc, :] > I2c[profc, :].max() / 2]
fy_width = x1[I2c[:, profc] > I2c[:, profc].max() / 2]
print(f"[6.3] 焦线:x 向半高宽 ≈ {(fx_width[-1]-fx_width[0])*1e3:.2f} mm(聚焦),"
      f"y 向 ≈ {(fy_width[-1]-fy_width[0])*1e3:.2f} mm(仍是衍射展宽)")

fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.8))
axs[0].imshow(I2c, origin="lower", extent=[x1[0], x1[-1], x1[0], x1[-1]],
              cmap="gray", vmax=I2c.max() * 0.3)
axs[0].set_title("Ex 6.3: cylindrical focus, z=zf=2000 m")
axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)")
# 6.2 的结果图放右边
axs[1].imshow(I2, origin="lower", extent=[x1[0], x1[-1], x1[0], x1[-1]],
              cmap="gray", vmax=I2.max() * 0.3)
axs[1].set_title(r"Ex 6.2: tilt($5\times10^{-5}$,45°)+focus(zf=4000)")
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("y (m)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_ex63.png"), dpi=150)
plt.close(fig)

# ================= 习题 6.5:波带片 =================
L1p = 0.05; Mp = 500; dx1p = L1p / Mp
x1p, _ = coords(L1p, Mp)
lam = 0.5e-6; wp = 6.25e-3; fp = 10.0
X1p, Y1p = np.meshgrid(x1p, x1p)
ta1 = zone_plate(X1p, Y1p, lam, fp, wp, binary=False)    # 式 (6.29)
ta2 = zone_plate(X1p, Y1p, lam, fp, wp, binary=True)     # 式 (6.30)
# (c) 采样判据:板边缘 r=w 处 chirp 局部频率 = r/(λf) = 1250 cyc/m < Nyquist 5000
print(f"[6.5] 板边缘局部频率 w/(λf) = {wp/(lam*fp):.0f} cyc/m < 1/(2Δx) = "
      f"{1/(2*dx1p):.0f} cyc/m(式 6.29 的采样判据满足,4 倍裕量)")
# (d) z=10 m:λz/(LΔx) = 1 → 临界采样,TF 与 IR 等价,用 TF
peaks1 = []; peaks2 = []
fig, axs = plt.subplots(2, 4, figsize=(13.6, 6.6))
for j, zt in enumerate((5.0, 8.0, 10.0, 12.0)):
    I_a = np.abs(propTF(ta1, L1p, lam, zt)) ** 2
    I_b = np.abs(propTF(ta2, L1p, lam, zt)) ** 2
    if zt == 10.0:
        I_lens = np.abs(propTF(lens_ta_plane(X1p, Y1p, lam, fp, wp),
                               L1p, lam, zt)) ** 2
        print(f"[6.5g] z=10 焦点峰值:cos 型 {I_a.max():.1f} / 二值型 {I_b.max():.1f}"
              f" / 理想透镜(同 f、同口径) {I_lens.max():.1f}")
        print(f"[6.5h] 透镜焦斑半高宽 "
              f"{fwhm(I_lens[Mp//2,:], x1p)*1e3:.3f} mm vs cos 型 "
              f"{fwhm(I_a[Mp//2,:], x1p)*1e3:.3f} mm —— 主瓣几乎同宽"
              f"(一级焦距相同、NA 相同),但波带片峰值低一个多量级,"
              f"且残留零阶与高阶焦点背景")
    for i, (Im, nm) in enumerate(((I_a, "cos plate"), (I_b, "binary plate"))):
        ax = axs[i, j]
        ax.imshow(Im ** 0.25, origin="lower",
                  extent=[x1p[0], x1p[-1], x1p[0], x1p[-1]], cmap="gray")
        ax.set_title(f"{nm}, z={zt:.0f} m", fontsize=9)
    peaks1.append(I_a.max()); peaks2.append(I_b.max())
print(f"[6.5e] 各距离峰值(cos/binary):"
      + ", ".join(f"z={zt:.0f}:{a:.1f}/{b:.1f}" for zt, a, b in
                  zip((5, 8, 10, 12), peaks1, peaks2)))
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_ex65.png"), dpi=150)
plt.close(fig)

# ================= 习题 6.7(a):相位光栅,零阶功率 =================
# 用无限周期模型(D1 = L1、L1/P = 100 整数,书 6.4.4 节):衍射阶退化为
# 样本 delta,零阶样本功率 = J0(m)^2(有限孔径模型的阶瓣跨多个样本,
# 峰值样本不等于整阶功率,不能用同样的口径比)。
D1g = 1e-2; P = 1e-4; L1g = 1e-2; Mg = 500; fg = 0.5
x1g, _ = coords(L1g, Mg)
X1g, Y1g = np.meshgrid(x1g, x1g)
rect2d = np.ones((Mg, Mg))                     # 光栅铺满阵列
print("[6.7] 相位光栅 exp(jm·cos(2πx/P)) 零阶功率分数(无限周期模型,"
      "理论 J0(m)^2,相位光栅总透射 = 1):")
for mm in (np.pi, np.pi / 2, np.pi / 4):
    u1g = np.exp(1j * mm * np.cos(2 * np.pi * X1g / P)) * rect2d
    u2g, L2g = propFF(u1g, L1g, lam, fg)
    dx2g = L2g / Mg
    Ig = np.abs(u2g) ** 2
    p0 = Ig[Mg // 2, Mg // 2] * dx2g ** 2 / L1g ** 2   # 零阶样本功率分数
    tot = Ig.sum() * dx2g ** 2 / L1g ** 2              # 总透射(相位光栅应 = 1)
    print(f"    m={mm:.4f}:零阶 {p0:.4f}(理论 {j0(mm)**2:.4f}),总透射 {tot:.4f}(=1)")

# ================= 习题 6.8:三角波光栅(1D) =================
P8 = (1 / 500) * 1e-3                          # 500 lines/mm → P = 2 um
D8 = 0.1e-3; lam8 = 0.5e-6; f8 = 0.5
dx8 = P8 / 8                                   # 提示:每个三角 8 个样本
M8 = 4096                                      # 偶数;L1 = M*dx8 = 1.024 mm
L18 = M8 * dx8
x18, _ = coords(L18, M8)
ux8 = tri_wave_1d(x18, P8)
u18 = ux8 * np.where(np.abs(x18) <= D8 / 2, 1.0, 0.0)
lf8 = lam8 * f8
u28 = np.sqrt(1 / lf8) * np.fft.ifftshift(
    np.fft.fft(np.fft.fftshift(u18))) * dx8
L28 = lf8 / dx8
x28, _ = coords(L28, M8)
I28 = np.abs(u28) ** 2
# 解析:(D1/2)·Σ sinc²(n/2)·sinc(D1(x2−nλf/P)/λf) → 幅度;阶位 nλf/P
u2a8 = np.zeros(M8)
for n in range(-15, 16):
    u2a8 += sinc(n / 2) ** 2 * sinc(D8 / lf8 * (x28 - n * lf8 / P8))
u2a8 = D8 / 2 / np.sqrt(lf8) * u2a8
I28a = np.abs(u2a8) ** 2
rel8 = np.max(np.abs(I28 - I28a)) / I28a.max()
c8 = int(np.argmax(I28))
print(f"[6.8] 三角光栅:P={P8*1e6:.0f} um(500 lines/mm),dx=P/8,"
      f"阶位 ±λf/P = ±{lf8/P8*1e3:.3f} mm;数值 vs 解析最大相对差 {rel8*100:.3f}%")

fig, ax = plt.subplots(figsize=(8.4, 3.9))
ax.plot(x28 * 1e3, I28 / I28.max(), label="discrete (1D)")
ax.plot(x28 * 1e3, I28a / I28a.max(), "--", label=r"analytic $\propto\mathrm{sinc}^2(n/2)$")
ax.set_xlim(-420, 420)
ax.set_xlabel("x (mm)"); ax.set_ylabel("normalized Irradiance"); ax.legend()
ax.set_title("Ex 6.8: triangle-wave grating (500 lines/mm), Fraunhofer 1D")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_ex68.png"), dpi=150)
plt.close(fig)

# ---- (c) 无限周期模型:L1/P = 1.024mm/2um = 512(整数)→ 阶为样本 delta ----
u18p = ux8                                   # 去掉 rect 窗(铺满阵列)
u28p = np.sqrt(1 / lf8) * np.fft.ifftshift(
    np.fft.fft(np.fft.fftshift(u18p))) * dx8
L28p = lf8 / dx8
x28p, _ = coords(L28p, M8)
I28p = np.abs(u28p) ** 2
Fp8 = I28p * (L28p / M8) / L18               # 1D 分数功率:I*dx2/L1(单位入射)
i0 = int(np.argmin(np.abs(x28p)))
i1 = int(np.argmin(np.abs(x28p - lf8 / P8)))
print(f"[6.8c] 无限周期模型:零阶 {Fp8[i0]:.4f}(理论 |c0|^2=0.25),"
      f"±1 阶 {Fp8[i1]:.4f}(理论 |(1/2)sinc^2(1/2)|^2=0.0411),"
      f"总透射 {Fp8.sum():.4f}(= <t^2>)")
plt.close(fig)
print("done ch6_exercises")
