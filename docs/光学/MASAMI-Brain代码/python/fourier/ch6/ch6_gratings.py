# -*- coding: utf-8 -*-
"""书 6.4 节:光栅三个模型(cos / 方波 / 一维 / 无限周期)数值验证。

复现:图 6.10(grating_cos 结果)、图 6.11(双波长分离)、
      图 6.13(grating_sqr 剖面)、图 6.14(一维模型)、
      图 6.16(无限周期模型的分数功率)。
验证:
  (1) grating_cos 数值 vs 式 (6.24) 解析(辐照度);
  (2) 阶位置 x = n·λf/P(零阶 + ±1 一阶 = ±2.5 mm;0.6 um 时 ±3.0 mm);
  (3) grating_sqr 数值 vs 式 (6.27) 解析(n=-5..5);
  (4) 无限周期模型分数功率:零阶 25%、±1 阶各 6.25%(书 6.4.4 节),
      并核对"入射功率分数 vs 透射功率分数"的口径差(笔记勘误);
  (5) 习题 6.6:M=1002(每周期 10.02 样本) vs 1100(整) 的落位影响。
输出:out_fig6_10/11/13/14/16.png。
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
from base import coords, sinc                  # noqa: E402  base.sinc = sin(πx)/(πx)
from ch5_prop import propFF                    # noqa: E402
from ch6_transmittance import square_wave_1d   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def run_grating_cos(lam=0.5e-6, M=500, D1=1.02e-3, P=1e-4, f=0.5, L1=1e-2):
    """书 6.4.1 节 grating_cos(式 6.20 + propFF + 式 6.24 解析)。"""
    dx1 = L1 / M
    x1, _ = coords(L1, M)
    X1, Y1 = np.meshgrid(x1, x1)
    u1 = 0.5 * (1 - np.cos(2 * np.pi * X1 / P)) * \
        np.where(np.abs(X1) <= D1 / 2, 1.0, 0.0) * \
        np.where(np.abs(Y1) <= D1 / 2, 1.0, 0.0)       # 式 (6.20)
    u2, L2 = propFF(u1, L1, lam, f)
    dx2 = L2 / M
    x2, _ = coords(L2, M)
    I2 = np.abs(u2) ** 2
    # 式 (6.24)(丢前导复指数,只算辐照度;书代码逐行对应)
    lf = lam * f
    X2, Y2 = np.meshgrid(x2, x2)
    u2a = (1 / lf) * D1 ** 2 / 2 * sinc(D1 / lf * Y2) * \
        (sinc(D1 / lf * X2) - 0.5 * sinc(D1 / lf * (X2 + lf / P))
         - 0.5 * sinc(D1 / lf * (X2 - lf / P)))
    I2a = np.abs(u2a) ** 2
    return x1, u1, x2, I2, I2a, lf, P, f, D1


# ---------- (1)(2) grating_cos ----------
x1, u1, x2, I2, I2a, lf, P, f, D1 = run_grating_cos()
prof = I2.shape[0] // 2
rel = np.max(np.abs(I2 - I2a)) / I2a.max()
order1 = lf / P
pk0 = I2[prof, prof]
pk1 = I2[prof, np.argmin(np.abs(x2 - order1))]
print(f"[1] grating_cos:数值 vs 式(6.24) 最大相对差 = {rel*100:.3f}%")
print(f"[2] 一阶位置:理论 ±λf/P = ±{order1*1e3:.2f} mm,实测 ±1 阶峰 "
      f"{I2[prof, np.argmin(np.abs(x2+order1))]:.2f} / 零阶 {pk0:.2f} / "
      f"{pk1:.2f}(书图 6.10(c):~1 / 4.2 / ~1)")

fig, axs = plt.subplots(1, 3, figsize=(12.4, 3.7))
axs[0].imshow(np.abs(u1) ** 2, origin="lower",
              extent=[x1[0], x1[-1], x1[0], x1[-1]], cmap="gray")
axs[0].set_title("source I1"); axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)")
axs[1].imshow(I2 ** (1 / 3), origin="lower",
              extent=[x2[0], x2[-1], x2[0], x2[-1]], cmap="gray")
axs[1].set_title("Fraunhofer I2^(1/3)"); axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("y (m)")
axs[2].plot(x2 * 1e3, I2[prof, :], label="discrete")
axs[2].plot(x2 * 1e3, I2a[prof, :], "--", label="analytic Eq.(6.24)")
axs[2].set_xlabel("x (mm)"); axs[2].set_ylabel("Irradiance"); axs[2].legend()
axs[2].set_title("grating_cos profile (cf. Fig. 6.10)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_10.png"), dpi=150)
plt.close(fig)

# ---------- 图 6.11:双波长分离 ----------
_, _, x2b, I2b, I2b_a, _, _, _, _ = run_grating_cos(lam=0.6e-6)
fig, ax = plt.subplots(figsize=(7.4, 3.9))
ax.plot(x2 * 1e3, I2[prof, :], label=r"0.5 $\mu$m")
ax.plot(x2b * 1e3, I2b[prof, :], "--", label=r"0.6 $\mu$m")
ax.set_xlabel("x (mm)"); ax.set_ylabel("Irradiance"); ax.legend()
ax.set_title("grating_cos, two wavelengths (cf. Fig. 6.11)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_11.png"), dpi=150)
plt.close(fig)
L2_05 = 0.5e-6 * f / (1e-2 / 500)              # 0.5um:propFF 观察面边长
L2_06 = 0.6e-6 * f / (1e-2 / 500)              # 0.6um:观察面边长(不同!)
print(f"[2b] 0.6um 一阶峰位 ±{0.6e-6*f/P*1e3:.2f} mm vs 0.5um ±{order1*1e3:.2f} mm"
      f"(两曲线 x2 网格不同:L2 = {L2_05:.4f} m vs {L2_06:.4f} m,不能用同一组 x2)")

# ---------- (3) grating_sqr(方波,式 6.25 + 式 6.27) ----------
lam = 0.5e-6
f = 0.5
P = 1e-4
D1 = 1e-3
L1 = 1e-2
M = 1000
dx1 = L1 / M
x1, _ = coords(L1, M)
print(f"[3] P/dx1 = {P/dx1:.0f}(整数),P/(2dx1) = {P/(2*dx1):.0f}(奇数)")
ux = square_wave_1d(x1, P)                     # rect ⊛ comb(书 1D 卷积)
X1, Y1 = np.meshgrid(x1, x1)
u1s = np.repeat(ux[None, :], M, axis=0)        # repmat(ux, M, 1)
u1s = u1s * np.where(np.abs(X1) <= D1 / 2, 1.0, 0.0) * \
    np.where(np.abs(Y1) <= D1 / 2, 1.0, 0.0)
u2s, L2s = propFF(u1s, L1, lam, f)
x2s, _ = coords(L2s, M)
I2s = np.abs(u2s) ** 2
prof_s = M // 2                                # y=0 行(本节 M=1000,与上节 prof 不同)
# 式 (6.27) 的 y=0 切片(n=-5..5,书代码逐行对应)
lfs = lam * f
u2a = np.zeros(M)
for n in range(-5, 6):
    u2a = u2a + sinc(n / 2) * sinc(D1 / lfs * (x2s - n * lfs / P))
u2a = D1 ** 2 / (2 * lfs) * u2a
I2s_a = np.abs(u2a) ** 2
rel_s = np.max(np.abs(I2s[prof_s, :] - I2s_a)) / I2s_a.max()
print(f"[3b] grating_sqr:数值 vs 式(6.27) 最大相对差 = {rel_s*100:.3f}%;"
      f"阶位 n·λf/P:±{lfs/P*1e3:.2f} mm(±1),±{3*lfs/P*1e3:.2f} mm(±3)")

fig, ax = plt.subplots(figsize=(7.4, 3.9))
ax.plot(x2s * 1e3, I2s[prof_s, :], label="discrete")
ax.plot(x2s * 1e3, I2s_a, "--", label="analytic Eq.(6.27)")
ax.set_xlabel("x (mm)"); ax.set_ylabel("Irradiance"); ax.legend()
ax.set_title("grating_sqr profile (cf. Fig. 6.13)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_13.png"), dpi=150)
plt.close(fig)

# ---------- grating_sqr1D(M=2200,书 6.4.3 节) ----------
M1 = 2200
dx1_1 = L1 / M1
x1_1, _ = coords(L1, M1)
print(f"[3c] 1D 模型:M={M1},P/dx1 = {P/dx1_1:.0f},L2 = λf/dx1 = "
      f"{lfs/dx1_1:.4f} m(2D M=1000 时 {lfs/dx1:.4f} m)")
ux1 = square_wave_1d(x1_1, P)
u1_1d = ux1 * np.where(np.abs(x1_1) <= D1 / 2, 1.0, 0.0)   # 只乘 rect(x1/D1)
u2_1d = np.sqrt(1 / lfs) * np.fft.ifftshift(
    np.fft.fft(np.fft.fftshift(u1_1d))) * dx1_1             # 书的三行(1D)
L2_1d = lfs / dx1_1
x2_1d, _ = coords(L2_1d, M1)
I2_1d = np.abs(u2_1d) ** 2

fig, axs = plt.subplots(2, 1, figsize=(8.2, 5.6))
axs[0].plot(x2_1d * 1e3, I2_1d / I2_1d.max())
axs[0].set_title("1D model, full vector (cf. Fig. 6.14a)")
axs[1].plot(x2_1d * 1e3, I2_1d / I2_1d.max())
axs[1].set_xlim(-10, 10)
axs[1].set_title("1D model, expanded (cf. Fig. 6.14b; 书正文误写 Fig. 6.11(b))")
for a in axs:
    a.set_ylabel("normalized Irradiance")
axs[1].set_xlabel("x (mm)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_14.png"), dpi=150)
plt.close(fig)

# ---------- (4) 无限周期模型(书 6.4.4 节,式 6.28) ----------
D1p = L1                                        # D1 = L1(铺满阵列)
L1p = 1e-2
Mp = 500
dx1p = L1p / Mp
x1p, _ = coords(L1p, Mp)
X1p, Y1p = np.meshgrid(x1p, x1p)
print(f"[4] L1/P = {L1p/P:.0f}(式 6.28 整数周期)")
u1p = 0.5 * (1 - np.cos(2 * np.pi * X1p / P)) * \
    np.where(np.abs(X1p) <= D1p / 2, 1.0, 0.0) * \
    np.where(np.abs(Y1p) <= D1p / 2, 1.0, 0.0)
u2p, L2p = propFF(u1p, L1p, lam, f)
dx2p = L2p / Mp
x2p, _ = coords(L2p, Mp)
I2p = np.abs(u2p) ** 2
Fpwr = I2p * dx2p ** 2 / (L1p ** 2)             # 书原代码
profp = Mp // 2
i0 = np.argmin(np.abs(x2p)); i1 = np.argmin(np.abs(x2p - lf / P))
print(f"[4b] 零阶 Fpwr = {Fpwr[profp, i0]:.4f}(书 25%),±1 阶 = "
      f"{Fpwr[profp, i1]:.4f}(书 6.25%)")
trans = (np.abs(u1p) ** 2).sum() * dx1p ** 2 / L1p ** 2
print(f"[4c] 振幅光栅真实透射功率 = {trans:.4f}(单位入射):"
      f"若除以透射功率,零阶 = {Fpwr[profp,i0]/trans*100:.1f}%、±1 阶各 "
      f"{Fpwr[profp,i1]/trans*100:.1f}% —— 书正文口径(除以透射功率)与代码矛盾")

fig, ax = plt.subplots(figsize=(7.0, 3.9))
ax.plot(x2p * 1e3, Fpwr[profp, :])
ax.set_xlabel("x (mm)"); ax.set_ylabel("Fractional Power")
ax.set_title("periodic model: fractional power (cf. Fig. 6.16)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "out_fig6_16.png"), dpi=150)
plt.close(fig)

# ---------- (5) 习题 6.6:M=1002 vs 1100 ----------
for Mt in (1000, 1002, 1100):
    dxt = L1 / Mt
    x1t, _ = coords(L1, Mt)
    uxt = square_wave_1d(x1t, P)
    X1t, Y1t = np.meshgrid(x1t, x1t)
    u1t = np.repeat(uxt[None, :], Mt, axis=0) * \
        np.where(np.abs(X1t) <= D1 / 2, 1.0, 0.0) * \
        np.where(np.abs(Y1t) <= D1 / 2, 1.0, 0.0)
    u2t, L2t = propFF(u1t, L1, lam, f)
    d = P / dxt
    fpk = (np.abs(u2t) ** 2).max() * (L2t / Mt) ** 2 / L1 ** 2
    print(f"[5] M={Mt}:每周期样本 P/dx1 = {d:.3f},最强阶样本的功率分数 = {fpk:.3e}"
          f"(M=1000 应为 2.5e-05;非整周期时阶模糊、峰值下降)")
print("done ch6_gratings")
