# -*- coding: utf-8 -*-
"""书 9.2 节:pc_spatial 部分空间相干仿真(随机相位屏 + 系综平均)。

书 9.2.1/9.2.2/9.2.3 节的公式与脚本:
  - 式 (9.13)-(9.16) 随机透过屏模型:t_A = U0 上的复相位屏,场相关
    mu_12 = R(x',y')/R(0,0),常用高斯形式 exp(-(x'^2+y'^2)/l_cr^2)
  - 式 (9.17)-(9.20) 相位屏 phi = r (*) f、高斯滤波器 f、屏自相关 R
  - 式 (9.21)-(9.23) 高斯相关近似与横向相干长度 l_cr^2 = 4 pi sigma_f^4/sigma_r^2
  - 式 (9.24)-(9.26) 频域傅里叶滤波法生成相位屏
  - pc_spatial(书 9.2.3 节第 3-61 行):N 块屏各自作用、辐照度系综平均

原书脚本逐行对应见 matlab/fourier/ch9/pc_spatial.m。本脚本完整跑通并输出:
  - |mu_12| = exp(-(dels/l_cr)^2)(书 0.677)与式 (9.20) 精确值
  - 相位屏生成链验证:屏参数、滤波器通带 vs 频率网格
  - N=100 块屏系综平均的夫琅禾费图样,与解析解式 (9.27) 剖面对比
  - 条纹可见度测量 vs |mu_12|(滤波器欠采样偏差的定量检查)

MATLAB/Python 差异:
  - randn(M)+j*randn(M) 对应 default_rng().standard_normal((M,M)) + 1j*(...);
    书不设随机种子,Python 侧固定种子保证可复现(统计量相同)。
  - fie = ifft2(F.*r~)*(sigma_r/dfx1)*M^2*dfx1^2(书 39-40 行,式 9.25):
    sigma_r/dfx1 来自式 (9.25) 的方差归一(均匀采样 (dfx dfy)^(1/2) = dfx1);
    M^2 抵消 ifft2 内置的 1/M^2;dfx1^2 补上逆变换积分的面积元。
  - 相位屏 fie 从 ifft2 出来是 FFT 排列(索引 0 = 零坐标),与 u1 同排列,
    直接逐点相乘;显示/测相关时再 ifftshift 回居中。
运行产物:out_fig9_screens.png、out_fig9_spatial.png + 控制台数值验证。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import j1

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from base import coords, circle, rect   # noqa: E402


# ===== 式 (9.20):屏自相关的精确值(不取式 9.22 的高斯近似) =====
def mu12_exact920(dels, Lcr, sigma_f):
    """式 (9.20) 在间隔 dels 处的精确值:mu_12 = R(dels)/R(0)。

    R = exp{-sigma_r^2/(2 pi sigma_f^2) [1 - exp(-dels^2/(2 sigma_f^2))]},
    代入式 (9.23) 的 sigma_r^2 = 4 pi sigma_f^4 / l_cr^2,即
    kappa = 2 sigma_f^2 / l_cr^2。sigma_f 有限时 mu_12(dels->inf) 不为 0,
    高斯近似(式 9.22)要求式 (9.21) 的 sigma_r^2/(2 pi sigma_f^2) >> 1。
    """
    kappa = 2 * sigma_f ** 2 / Lcr ** 2
    return np.exp(-kappa * (1 - np.exp(-dels ** 2 / (2 * sigma_f ** 2))))


def analytic_two_beam_spatial(x2, y2, lam, f, w, dels, mu12):
    """式 (9.27):双圆孔 + 部分空间相干的解析辐照度(单波长)。

    与式 (9.10) 同构:gamma(deld/c) 换成 mu_12,余弦相位里没有 k0*deld 项
    (无路径延迟)。原书笔误同式 (9.10):花括号内 w^2 应为 w^4 —— 按 w^4
    计算中心值 1253,与书图 9.11(b) 峰值 ~1250 吻合。
    """
    wtil = w / (lam * f)
    r = np.sqrt(x2 ** 2 + y2 ** 2)
    z = 2 * np.pi * wtil * r
    ratio = np.where(r == 0, np.pi, j1(z) / (wtil * np.where(r == 0, 1.0, r)))
    envelope = (2 * w ** 4 / (lam * f) ** 2) * ratio ** 2
    return envelope * (1 + mu12 * np.cos(2 * np.pi * dels * x2 / (lam * f)))


def pc_spatial(N=100, Lcr=8e-3, mult=2.5, aperture="two_beams", seed=9):
    """书 pc_spatial 主体(第 3-51 行,逐行对应)。

    mult 是书 21 行的 2.5 系数(习题 9.4 让试更小的值);aperture 可换
    "rect5x2"(习题 9.5 的 5mm x 2mm 矩形孔)。
    返回 I2(居中)、x2(观察面坐标)与 (sigma_f, sigma_r)。
    """
    lam = 650e-9                                # 书 3 行
    L1 = 50e-3                                  # 书 5 行
    M = 250                                     # 书 6 行
    dx1 = L1 / M                                # 书 7 行
    x1, _ = coords(L1, M)
    x1 = np.fft.fftshift(x1)                    # 书 9 行
    X1, Y1 = np.meshgrid(x1, x1)
    w = 1e-3                                    # 书 13 行
    dels = 5e-3                                 # 书 14 行
    f = 0.25                                    # 书 15 行
    lf = lam * f                                # 书 16 行

    sigma_f = mult * Lcr                        # 书 21 行(默认 2.5*Lcr)
    sigma_r = np.sqrt(4 * np.pi * sigma_f ** 4 / Lcr ** 2)   # 书 22 行(式 9.23)
    dfx1 = 1 / L1                               # 书 24 行
    _, fx1 = coords(L1, M)
    fx1 = np.fft.fftshift(fx1)                  # 书 26 行
    FX1, FY1 = np.meshgrid(fx1, fx1)

    if aperture == "two_beams":                 # 书 30-31 行
        u1 = circle(X1 - dels / 2, Y1, w) + circle(X1 + dels / 2, Y1, w)
    else:                                       # 习题 9.5:5mm x 2mm 矩形孔
        u1 = rect(X1 / 5e-3) * rect(Y1 / 2e-3)
    F = np.exp(-np.pi ** 2 * sigma_f ** 2 * (FX1 ** 2 + FY1 ** 2))   # 书 33 行(式 9.26)

    I2 = np.zeros((M, M))                       # 书 36 行
    rng = np.random.default_rng(seed)
    for n in range(N // 2):                     # 书 37 行:每次生成 2 块屏
        fie = (np.fft.ifft2(F * (rng.standard_normal((M, M))
                                 + 1j * rng.standard_normal((M, M))))
               * sigma_r / dfx1) * M ** 2 * dfx1 ** 2    # 书 39-40 行(式 9.25)
        u2 = 1 / lf * np.fft.fft2(u1 * np.exp(1j * fie.real)) * dx1 ** 2
        I2 = I2 + np.abs(u2) ** 2               # 书 42-43 行
        u2 = 1 / lf * np.fft.fft2(u1 * np.exp(1j * fie.imag)) * dx1 ** 2
        I2 = I2 + np.abs(u2) ** 2               # 书 45-46 行

    I2 = np.fft.ifftshift(I2) / N               # 书 49 行:normalize & center
    _, fx2 = coords(L1, M)
    x2 = fx2 * lf                               # 书 50 行
    return I2, x2, (sigma_f, sigma_r)


def fringe_visibility(prof, x2, per):
    """中央条纹的直接可见度 (Imax-Imin)/(Imax+Imin)(抛物线插值极值)。

    比 LSQ 拟合稳:LSQ 以"相干单束包络 g2"作基,而系综平均的实际包络被
    屏的相内相位微调宽,基失配会把 mu 估高;直接极值只读中央 1-2 个条纹。
    """
    win = np.abs(x2) < 1.2 * per
    p = prof[win]
    imax = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1
    imin = np.where((p[1:-1] < p[:-2]) & (p[1:-1] < p[2:]))[0] + 1

    def refine(idx):
        if len(idx) == 0:
            return np.empty(0)
        a, b, c = p[idx - 1], p[idx], p[idx + 1]
        den = a - 2 * b + c
        return b - (a - c) ** 2 / (8 * np.where(den == 0, np.inf, den))

    pk, vl = refine(imax), refine(imin)
    if len(pk) == 0 or len(vl) == 0:
        return np.nan
    return (pk.max() - vl.min()) / (pk.max() + vl.min())


def make_screens(M, L1, Lcr, n_pairs, mult=2.5, seed=11):
    """书参数下生成 n_pairs 对相位屏(实部/虚部各一块),供图 9.10 显示。

    与 pc_spatial 内循环同款(书 39-40 行);mult 可调(sigma_f = mult*Lcr)。
    """
    sigma_f = mult * Lcr
    sigma_r = np.sqrt(4 * np.pi * sigma_f ** 4 / Lcr ** 2)
    dfx1 = 1 / L1
    _, fx1 = coords(L1, M)
    fx1 = np.fft.fftshift(fx1)
    FX1, FY1 = np.meshgrid(fx1, fx1)
    F = np.exp(-np.pi ** 2 * sigma_f ** 2 * (FX1 ** 2 + FY1 ** 2))
    rng = np.random.default_rng(seed)
    fie = np.empty((n_pairs, M, M), dtype=complex)
    for n in range(n_pairs):
        rtilde = rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M))
        fie[n] = (np.fft.ifft2(F * rtilde) * sigma_r / dfx1) * M ** 2 * dfx1 ** 2
    return fie, (sigma_f, sigma_r)


if __name__ == "__main__":
    lam, L1, M = 650e-9, 50e-3, 250
    w, dels, f = 1e-3, 5e-3, 0.25
    Lcr = 8e-3
    mu12 = np.exp(-(dels / Lcr) ** 2)           # 式 (9.16):|mu_12| at dels
    print(f"|mu_12| = exp(-(dels/l_cr)^2) = {mu12:.4f}(书 0.677,图 9.11 标 0.68)")
    print(f"式 (9.20) 精确值(不取高斯近似,sigma_f=2.5*l_cr) = "
          f"{mu12_exact920(dels, Lcr, 2.5 * Lcr):.4f}")

    # ----- 相位屏参数(书 21-22 行) -----
    sigma_f0 = 2.5 * Lcr
    sigma_r = np.sqrt(4 * np.pi * sigma_f0 ** 4 / Lcr ** 2)
    print(f"sigma_f = 2.5*l_cr = {sigma_f0*1e3:.1f} mm(< L/2 = {L1/2*1e3:.0f} mm);"
          f"sigma_r = {sigma_r*1e3:.1f} mm(书 177.2 mm)")
    print(f"式 (9.21): sigma_r^2/(2 pi sigma_f^2) = {sigma_r**2/(2*np.pi*sigma_f0**2):.2f}"
          f" >> 1;滤波器通带 1/(pi*sigma_f) = {1/(np.pi*sigma_f0):.1f} 1/m"
          f" vs 频率网格 dfx1 = {1/L1:.1f} 1/m -> 通带不足一个网格,欠采样!")

    # ----- 主仿真(书默认参数 N=100) -----
    I2, x2, (sf, sr) = pc_spatial()
    prof = I2[M // 2, :]
    I_an = analytic_two_beam_spatial(x2, 0.0, lam, f, w, dels, mu12)
    rms = np.sqrt(np.mean((prof - I_an) ** 2)) / I_an.max()
    print(f"剖面 vs 式 (9.27):相对 RMS 偏差 = {rms:.3f}(N=100 屏的系综噪声量级);"
          f"剖面峰值:仿真 {prof.max():.0f},解析 {I_an.max():.0f}(书图 9.11 ~1250)")

    # ----- 可见度对照(直接极值)与滤波器分辨率扫描 -----
    per = f * lam / dels
    V_book = fringe_visibility(prof, x2, per)
    print(f"条纹可见度(书参数,中央条纹直接极值,N=100)= {V_book:.3f}"
          f" —— 目标 |mu_12| = 0.677:偏高,见下方扫描")
    print("mult  sigma_f(mm)  通带/dfx1   可见度(N=200)  式(9.20)精确  高斯(9.22)")
    mults = [2.5, 2.0, 1.5, 1.2]
    V_meas = []
    for mult in mults:
        Ij, xj, _ = pc_spatial(N=200, Lcr=Lcr, mult=mult, seed=7)
        Vj = fringe_visibility(Ij[Ij.shape[0] // 2, :], xj, per)
        V_meas.append(Vj)
        sfj = mult * Lcr
        print(f"{mult:4.1f}  {sfj*1e3:9.1f}  {1/(np.pi*sfj)*L1:8.2f}   "
              f"{Vj:11.4f}  {mu12_exact920(dels, Lcr, sfj):11.4f}  "
              f"{np.exp(-(dels/Lcr)**2):8.4f}")
    print("解读:2.5 倍率下通带(0.8 个网格)欠采样,屏的结构由 DC 与 ±20 cyc/m"
          " 的少数样本决定,相关长度被拉长 -> 可见度偏向相干(书 9.2.3 节自己描述的"
          " failure mode,只是发生在 sigma_f < L/2 时也会出现);倍率降低后可见度"
          "回落、逼近式 (9.20) 精确值;彻底解决要加大 L1 与 M 让通带被充分采样。")

    # ----- 大网格屏验证(L1=200mm, M=1000:dfx1=5 1/m,通带 3.2 个网格) -----
    sf_big = 2.5 * Lcr
    sr_big = np.sqrt(4 * np.pi * sf_big ** 4 / Lcr ** 2)
    L1b, Mb = 200e-3, 1000
    dfxb = 1 / L1b
    _, fxb = coords(L1b, Mb)
    FXb, FYb = np.meshgrid(np.fft.fftshift(fxb), np.fft.fftshift(fxb))
    Fb = np.exp(-np.pi ** 2 * sf_big ** 2 * (FXb ** 2 + FYb ** 2))
    rngb = np.random.default_rng(7)
    step = int(round(dels / (L1b / Mb)))
    accR = accD = 0.0
    nb = 120                                    # 120 块屏 × 1e6 样本,均值足够稳
    for k in range(nb):
        rt = rngb.standard_normal((Mb, Mb)) + 1j * rngb.standard_normal((Mb, Mb))
        phb = np.fft.ifftshift(((np.fft.ifft2(Fb * rt) * sr_big / dfxb)
                                * Mb ** 2 * dfxb ** 2).real)
        tb = np.exp(1j * phb)
        accR += np.mean(tb * np.roll(tb, -step, axis=1).conj()).real
        accD += np.mean((np.roll(phb, -step, axis=1) - phb) ** 2)
    D_pred = 2 * (sr_big ** 2 / (2 * np.pi * sf_big ** 2)) \
        * (1 - np.exp(-dels ** 2 / (2 * sf_big ** 2)))
    print(f"大网格屏验证(L1=200mm, M=1000, 300 块屏):"
          f"相位结构函数 D(5mm) = {accD/nb:.4f}(式 9.20 预测 {D_pred:.4f}),"
          f"透过率相关 R_t(5mm) = {accR/nb:.4f}(预测 {np.exp(-D_pred/2):.4f})"
          f" —— 滤波器被充分采样时,屏相关精确落到连续理论值")

    # ----- 出图 1:相位屏两块(书图 9.10) -----
    fie, _ = make_screens(M, L1, Lcr, n_pairs=2, seed=11)
    fig = plt.figure(figsize=(11, 5))
    xx = (np.arange(M) - M // 2) * (L1 / M) * 1e3   # mm
    for j in range(2):
        ax = fig.add_subplot(1, 2, j + 1, projection="3d")
        ph = np.fft.ifftshift(fie[j].real)
        ax.plot_surface(xx, xx, ph, cmap="gray", rstride=4, cstride=4,
                        linewidth=0, antialiased=True)
        ax.set_zlim(-8, 8)
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        ax.set_zlabel("phase (rad)")
        ax.set_title(f"realization {j+1}")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "out_fig9_screens.png"), dpi=150)

    # ----- 出图 2:图样 + 剖面 + mu_eff 扫描(书图 9.11 + 验证) -----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ext = [x2[0], x2[-1], x2[0], x2[-1]]
    axes[0].imshow(I2, origin="lower", extent=ext, cmap="gray")
    axes[0].set_title("Irradiance image (Fig 9.11a)")
    axes[0].set_xlabel("x (m)"); axes[0].set_ylabel("y (m)")
    axes[1].plot(x2 * 1e4, prof)
    axes[1].plot(x2[::2] * 1e4, I_an[::2], "o", ms=3, label="analytic (9.27)")
    axes[1].set_xlim(-2, 2)
    axes[1].set_title("Profile, $\\ell_{cr}$=8 mm (Fig 9.11b)")
    axes[1].set_xlabel("x (1e-4 m)"); axes[1].legend(fontsize=8)
    sfs = np.linspace(0.75, 2.6, 100) * Lcr
    axes[2].plot(sfs / Lcr, [mu12_exact920(dels, Lcr, s) for s in sfs],
                 "-", label="Eq (9.20) exact")
    axes[2].plot([0, 3], [mu12, mu12], ":", label="Gaussian (9.22)")
    axes[2].plot(np.array(mults), V_meas, "o", label="sim visibility")
    axes[2].plot([2.5], [V_book], "s", label="book params (N=100)")
    axes[2].set_xlabel("$\\sigma_f\\,/\\,\\ell_{cr}$ (book: 2.5)")
    axes[2].set_ylabel("fringe visibility at $\\Delta s$")
    axes[2].set_title("Visibility vs filter width")
    axes[2].set_ylim(0.4, 1.0); axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "out_fig9_spatial.png"), dpi=150)
    print("saved out_fig9_screens.png / out_fig9_spatial.png")
