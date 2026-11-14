# -*- coding: utf-8 -*-
"""书 9.1 节:pc_temp 部分时间相干仿真(双光束杨氏干涉 + 高斯谱线型)。

书 9.1.1/9.1.2/9.1.3 节的公式与脚本:
  - 式 (9.1) 归一化高斯谱线型、式 (9.2) FWHM 宽度参数 b
  - 式 (9.3)/(9.4) 复时间度和 gamma(tau) = FT{S^(nu)}(正频率)
  - 式 (9.6)/(9.7) 相干时间 tau_c = 0.664/dnu、相干长度 l_c = c*tau_c
  - 式 (9.8)/(9.9) 总辐照度 = 谱加权辐照度的积分/离散和
  - pc_temp(书 9.1.3 节第 3-57 行):N 条谱线各自传播、按 S*dnu 加权累加

原书脚本逐行对应见 matlab/fourier/ch9/pc_temp.m。本脚本完整跑通并输出:
  - tau_c / l_c / |gamma(deld/c)| 解析值(书 3.32e-10 s、10 cm、0.6733)
  - 离散谱线的 gamma 数值和 vs 连续解析(式 9.3 的 DFT 近似)
  - N=51 条谱线叠加的夫琅禾费图样,与解析解式 (9.10) 剖面对比
  - deld = 0/5/10/50 cm 系列(书图 9.6,|gamma| = 1.0/0.67/0.21/0)

MATLAB/Python 差异:
  - 书第 18-19 行:x1 升序居中 -> fftshift(x1) 把"坐标数组本身"重排为 FFT
    排列,场数组 u 天生 DC 在索引 0,循环里的 fft2 不用再移位(省时);
    numpy 的 fftshift/ifftshift 行为与 MATLAB 一致(偶数 M 时两者相同)。
  - circ(sqrt(...)/w) 是半径 w 的圆盘(Goodman 半径约定,见 ch1 勘误),
    Python 用 base.circle(X1, Y1, w) 直接传半径。
  - 书第 56 行剖面 I2(M/2+1,:)(1-based)对应 I2[M//2, :](0-based)。
  - 书不设随机种子(本题无随机数);Python 侧固定种子仅为复现一致。
运行产物:out_fig9_temporal.png、out_fig9_temporal_series.png + 控制台验证。
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
from base import coords, circle, rect   # noqa: E402  公共积木(ch1 起)


# ===== 高斯谱线型与复时间度和(式 9.1-9.4,书 9.1.1 节) =====
def gaussian_lineshape(nu, nu0, b):
    """归一化高斯谱线型 S^(nu) = 1/(sqrt(pi) b) exp(-((nu-nu0)/b)^2)(式 9.1)。

    全频率积分为 1;FWHM 线宽 dnu 与宽度参数的关系是式 (9.2):
    b = dnu/(2 sqrt(ln2))(半高值 1/2 代入解出)。
    """
    return 1 / (np.sqrt(np.pi) * b) * np.exp(-((nu - nu0) / b) ** 2)


def gamma_temporal(tau, delnu, nu0):
    """复时间度和 gamma(tau)(式 9.4):
    gamma = exp[-(pi dnu tau/(2 sqrt(ln2)))^2] * exp(-j 2 pi nu0 tau)。

    高斯包络决定模长 |gamma|(式 9.11 的形式),载波相位由中心频率给出:
    gamma(tau) = FT{S^} 在正频率上的积分(式 9.3)。
    """
    return np.exp(-(np.pi * delnu * tau / (2 * np.sqrt(np.log(2)))) ** 2) \
        * np.exp(-1j * 2 * np.pi * nu0 * tau)


def analytic_two_beam(x2, y2, lam0, f, w, dels, deld, gamma_abs):
    """式 (9.10):双圆孔 + 延迟 deld 的解析辐照度(时域部分相干情形)。

    I2 = 2 w^4/(lam0 f)^2 * [J1(z)/(w r/(lam0 f))]^2
         * [1 + |gamma(deld/c)| cos(2 pi dels x2/(lam0 f) + k0 deld)]
    其中 z = 2 pi w r/(lam0 f),r = sqrt(x2^2+y2^2);r->0 时
    J1(z)/(w r/(lam0 f)) -> pi(用极限值避免除零)。

    原书笔误:式 (9.10) 花括号内印的是 w^2,应为 w^4 —— 按 w^2 计算中心
    辐照度是 7.5e8(比仿真大 6 个量级);按 w^4 得 ~1.2e3,与书图 9.5(b)
    峰值 ~1200 吻合(等价形式 2 (pi w^2/(lam0 f))^2 [2 J1(z)/z]^2)。
    """
    wtil = w / (lam0 * f)                       # w/(lam0 f)
    r = np.sqrt(x2 ** 2 + y2 ** 2)
    z = 2 * np.pi * wtil * r                    # J1 的宗量
    ratio = np.where(r == 0, np.pi, j1(z) / (wtil * np.where(r == 0, 1.0, r)))
    envelope = (2 * w ** 4 / (lam0 * f) ** 2) * ratio ** 2
    k0 = 2 * np.pi / lam0
    return envelope * (1 + gamma_abs * np.cos(2 * np.pi * dels * x2 / (lam0 * f)
                                              + k0 * deld))


def pc_temp(deld=5e-2, delnu=2e9, lineshape="gauss", amps=(1.0, 1.0), N=None):
    """书 pc_temp 主体(第 3-48 行,逐行对应)。

    lineshape="rect" 供习题 9.1 复用(式 9.29);amps=(A1,A2) 供习题 9.3
    复用(两束振幅不等,书取缺省 (1,1),行为不变)。
    返回 I2(居中排列)、x2(观察面坐标,居中升序)与 (nu, S, dnu) 谱线数组。
    """
    lam0 = 650e-9                               # 书 3 行:center wavelength
    c = 3e8                                     # 书 4 行
    nu0 = c / lam0                              # 书 6 行:center frequency
    N = 51 if N is None else N                  # 书 9 行:number of components (odd)
    b = delnu / (2 * np.sqrt(np.log(2)))        # 书 11 行:FWHM scaling(式 9.2)
    dnu = 4 * delnu / N                         # 书 12 行:freq interval
    L1 = 50e-3                                  # 书 15 行:source plane side length
    M = 250                                     # 书 16 行
    dx1 = L1 / M                                # 书 17 行
    x1, _ = coords(L1, M)                       # 书 18 行:居中升序
    x1 = np.fft.fftshift(x1)                    # 书 19 行:坐标移位 -> 场天生 FFT 排列
    X1, Y1 = np.meshgrid(x1, x1)                # 书 20 行
    w = 1e-3                                    # 书 23 行:radius
    dels = 5e-3                                 # 书 24 行:transverse separation
    f = 0.25                                    # 书 26 行:focal dist for Fraunhofer
    lf = lam0 * f                               # 书 27 行

    I2 = np.zeros((M, M))                       # 书 30 行
    nu_arr = np.empty(N)
    S_arr = np.empty(N)
    for n in range(1, N + 1):                   # 书 31 行
        nu = (n - (N + 1) / 2) * dnu + nu0      # 书 33 行
        if lineshape == "gauss":                # 书 34 行
            S = gaussian_lineshape(nu, nu0, b)
        else:                                   # 习题 9.1:矩形线型(式 9.29)
            S = rect((nu - nu0) / delnu) / delnu
        k = 2 * np.pi * nu / c                  # 书 35 行:逐谱线更新 k!
        # 书 37-39 行:两束圆盘,延迟项 exp(j k deld) 乘在 x=-dels/2 一束上
        # (习题 9.3:amps 给两束不等振幅 A1、A2)
        u = amps[0] * circle(X1 - dels / 2, Y1, w) \
            + amps[1] * circle(X1 + dels / 2, Y1, w) * np.exp(1j * k * deld)
        # 书 41 行:Fraunhofer 图样(忽略前端相位;中心波长定标)
        u2 = 1 / lf * np.fft.fft2(u) * dx1 ** 2
        I2 = I2 + S * np.abs(u2) ** 2 * dnu     # 书 43 行:式 (9.9) 的离散近似
        nu_arr[n - 1], S_arr[n - 1] = nu, S

    I2 = np.fft.ifftshift(I2)                   # 书 46 行:normalize/center
    _, fx2 = coords(L1, M)                      # 书 47 行:-1/(2dx1)..1/(2dx1)-1/L1
    x2 = fx2 * lf                               # 观察面坐标(中心波长定标)
    return I2, x2, (nu_arr, S_arr, dnu, lam0, f, w, dels)


def fringe_visibility(prof, x2, per):
    """中央条纹的直接可见度 (Imax-Imin)/(Imax+Imin),抛物线插值极值。

    等幅双束时 V 应逼近 |gamma|(式 9.31 的 A1=A2 情形);包络在一个条纹
    周期内略有衰减(空气瓣),直接极值会略微高估,故只在中心 ±1 周期内取。
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


if __name__ == "__main__":
    # ----- 参数与相干量(书 9.1.1/9.1.3 节数字) -----
    lam0, c = 650e-9, 3e8
    nu0 = c / lam0
    delnu = 2e9
    deld = 5e-2
    tau_c = 0.664 / delnu                       # 式 (9.6)
    l_c = c * tau_c                             # 式 (9.7)
    print(f"tau_c = 0.664/dnu = {tau_c:.3e} s(书 3.32e-10);"
          f"l_c = c*tau_c = {l_c*1e2:.1f} cm(书 10 cm)")
    tau = deld / c                              # tau = deld/c = 1.67e-10 s
    g_abs = abs(gamma_temporal(tau, delnu, nu0))
    print(f"|gamma(deld/c)| 解析 = {g_abs:.4f}(书 0.6733)")

    I2, x2, (nu_arr, S_arr, dnu, lam0_, f_, w_, dels_) = pc_temp(deld=deld)
    M = I2.shape[0]

    # ----- 数值验证 1:离散谱线的 gamma 数值和(式 9.3 的 DFT 近似) -----
    print(f"谱线权重和 sum(S*dnu) = {np.sum(S_arr) * dnu:.6f}(归一化线型应 = 1)")
    gam_num = np.sum(S_arr * np.exp(-2j * np.pi * nu_arr * tau)) * dnu
    print(f"离散度和 |sum S e^-j2pi nu tau dnu| = {abs(gam_num):.4f} "
          f"(解析 {g_abs:.4f},相对差 {abs(abs(gam_num) - g_abs) / g_abs:.2e})")

    # ----- 数值验证 2:剖面 vs 式 (9.10)(w^4 勘误版) -----
    I_an = analytic_two_beam(x2, 0.0, lam0_, f_, w_, dels_, deld, g_abs)
    prof = I2[M // 2, :]                        # 书 56 行:I2(M/2+1,:)
    rel = np.max(np.abs(prof - I_an)) / I_an.max()
    env_pk = 2 * np.pi ** 2 * w_ ** 4 / (lam0_ * f_) ** 2
    print(f"剖面峰值:仿真 {prof.max():.1f},解析 {I_an.max():.1f}"
          f"(书图 9.5(b) 峰值 ~1200;包络中心值 2 pi^2 w^4/(lam0 f)^2 = {env_pk:.0f},"
          f"乘 (1+|gamma|cos(k0 deld)) 的条纹峰 = {env_pk * (1 + g_abs):.0f})")
    print(f"剖面 vs 式 (9.10):最大相对偏差 = {rel:.2e}(峰值归一)"
          f" —— 像元化圆盘的边缘采样差异,即书图 9.6 评注的 lumpiness 量级")
    # 反证:若按书印的 w^2 计算,解析峰值会大 1/w^2 = 1e6 倍
    print(f"反证勘误:若按书印 w^2 计算,解析中心值 = {env_pk / w_ ** 2:.2e}"
          f"(比仿真大 6 个量级,式 (9.10) 的 w^2 应为 w^4)")

    # ----- 数值验证 3:条纹可见度 = |gamma| -----
    per = f_ * lam0_ / dels_                    # 条纹周期 lam0 f / dels
    V = fringe_visibility(prof, x2, per)
    print(f"条纹可见度(中央 ±1 周期,插值极值)= {V:.4f}(目标 |gamma| = {g_abs:.4f})")

    # ----- 出图 1:书图 9.5 + 图 9.2 + 图 9.3 复现 -----
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.6))
    ext = [x2[0], x2[-1], x2[0], x2[-1]]
    axes[0].imshow(I2, origin="lower", extent=ext, cmap="gray")
    axes[0].set_title("Irradiance image (Fig 9.5a)")
    axes[0].set_xlabel("x (m)"); axes[0].set_ylabel("y (m)")
    axes[1].plot(x2, prof)
    axes[1].plot(x2[::3], I_an[::3], "o", ms=3, label="analytic (9.10)")
    axes[1].set_xlim(-2e-4, 2e-4); axes[1].set_ylim(0, 1500)
    axes[1].set_title("Profile, $\\Delta d$=5 cm (Fig 9.5b)")
    axes[1].set_xlabel("x (m)"); axes[1].legend(fontsize=8)
    tq = np.linspace(-2 / delnu, 2 / delnu, 400)
    axes[2].plot(tq * 1e9, np.abs(np.exp(-(np.pi * delnu * tq / (2 * np.sqrt(np.log(2)))) ** 2)))
    axes[2].plot([tau * 1e9], [g_abs], "o")
    axes[2].set_title("|gamma(tau)|, tau=0.167 ns (Fig 9.2)")
    axes[2].set_xlabel("tau (ns)"); axes[2].set_ylim(0, 1)
    axes[3].plot((nu_arr - nu0) / 1e9, S_arr * dnu, "o-", ms=3)
    axes[3].set_title("Sampled PSD, 51 lines (Fig 9.3)")
    axes[3].set_xlabel("(nu-nu0)/dnu (GHz)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "out_fig9_temporal.png"), dpi=150)

    # ----- 出图 2:deld = 0/5/10/50 cm 系列(书图 9.6) -----
    fig, axes = plt.subplots(2, 4, figsize=(15, 6))
    for j, dd in enumerate([0.0, 5e-2, 1e-1, 5e-1]):
        Ij, xj, _ = pc_temp(deld=dd)
        g_j = np.exp(-(np.pi * delnu * (dd / c) / (2 * np.sqrt(np.log(2)))) ** 2)
        ext = [xj[0], xj[-1], xj[0], xj[-1]]
        axes[0, j].imshow(Ij, origin="lower", extent=ext, cmap="gray")
        axes[0, j].set_title(f"$\\Delta d$={dd*1e2:.0f} cm, |gamma|={g_j:.2f}")
        axes[1, j].plot(xj * 1e6, Ij[Ij.shape[0] // 2, :])
        axes[1, j].set_xlim(-200, 200)          # 书图 9.6:放大到中央 ±2e-4 m
        axes[1, j].set_xlabel("x (um)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "out_fig9_temporal_series.png"), dpi=150)
    print("saved out_fig9_temporal.png / out_fig9_temporal_series.png")
