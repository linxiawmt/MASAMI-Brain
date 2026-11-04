# -*- coding: utf-8 -*-
"""书 9.1 节:pc_temp 部分时间相干仿真(双光束杨氏干涉 + 高斯谱线型)。

原书脚本逐行对应见 matlab/fourier/ch9/pc_temp.m;本脚本完整跑通并输出:
  - 相干时间/相干长度(式 9.6/9.7)与 |gamma(tau)| 解析值(式 9.4/9.11)
  - 离散谱线的 gamma 数值和(S dnu 权重和、复数度和的模)与解析对比
  - N 条谱线叠加的夫琅禾费图样(式 9.9 的离散实现)
  - 与解析解式 (9.10) 的剖面逐点对比 + 干涉条纹可见度测量(= |gamma|)

MATLAB/Python 差异:
  - 书第 18-19 行:x1 升序居中 -> fftshift(x1) 把"坐标数组本身"重排为 FFT
    排列,场数组 u 天生 DC 在索引 0(省掉循环内的 fftshift);numpy 行为一致。
  - 书第 46 行 ifftshift(I2) 把累加的辐照度移回居中;numpy 相同。
  - circ(sqrt(...)/w) 是半径 w 的圆盘(Goodman 约定,见 ch1 勘误),Python 用
    base.circle(X1, Y1, w) 直接传半径。
  - 书第 56 行剖面 I2(M/2+1,:) (1-based) 对应 I2[M//2, :](0-based)。
运行产物:out_fig9_temporal.png、out_fig9_temporal_series.png + 控制台数值验证。
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
sys.path.insert(0, HERE)            # python -I 隔离模式不含脚本目录,显式加入
from base import coords, circle, rect   # noqa: E402


# ===== 高斯谱线型与复时间度和(式 9.1-9.4,书 9.1.1 节) =====
def gaussian_lineshape(nu, nu0, b):
    """归一化高斯谱线型 S^(nu) = 1/(sqrt(pi)b) exp(-((nu-nu0)/b)^2)(式 9.1)。"""
    return 1 / (np.sqrt(np.pi) * b) * np.exp(-((nu - nu0) / b) ** 2)


def gamma_temporal(tau, delnu, nu0):
    """复时间度和 gamma(tau)(式 9.4):
    gamma = exp[-(pi dnu tau/(2 sqrt(ln2)))^2] * exp(-j 2 pi nu0 tau)。
    高斯包络决定模长 |gamma|(式 9.11),载波相位由中心频率给出;
    delnu = FWHM(式 9.2:b = dnu/(2 sqrt(ln2)))。
    """
    return np.exp(-(np.pi * delnu * tau / (2 * np.sqrt(np.log(2)))) ** 2) \
        * np.exp(-1j * 2 * np.pi * nu0 * tau)


def analytic_two_beam(x2, y2, lam0, f, w, dels, deld, c, gamma_abs):
    """式 (9.10):双圆孔 + 延迟 deld 的解析辐照度(时域部分相干情形)。

    I2 = 2/(lam0 f)^2 * [w^2 J1(2pi w r/(lam0 f)) / (w r/(lam0 f))]^2
         * [1 + |gamma(deld/c)| cos(2pi dels x2/(lam0 f) + k0 deld)]
    r -> 0 时 J1(z)/w~r -> pi(z = 2 pi w r/(lam0 f)),用极限值避免除零。
    """
    wtil = w / (lam0 * f)                       # w/(lam0 f)
    r = np.sqrt(x2 ** 2 + y2 ** 2)
    z = 2 * np.pi * wtil * r                    # J1 的宗量
    ratio = np.where(r == 0, np.pi, j1(z) / (wtil * np.where(r == 0, 1.0, r)))
    envelope = (2 / (lam0 * f) ** 2) * (w ** 2 * ratio) ** 2
    k0 = 2 * np.pi / lam0
    return envelope * (1 + gamma_abs * np.cos(2 * np.pi * dels * x2 / (lam0 * f)
                                              + k0 * deld))


def pc_temp(deld=5e-2, delnu=2e9, lineshape="gauss", N=None, amps=(1.0, 1.0)):
    """书 pc_temp 主体(第 3-48 行)。

    lineshape="rect" 供习题 9.1 复用(式 9.29);amps=(A1,A2) 供习题 9.3
    复用(两束振幅不等,书代码取缺省 (1,1),行为不变)。
    返回 I2(居中排列)、x2(观察面坐标,居中升序)与谱线数组(nu, S)。
    """
    lam0 = 650e-9
    c = 3e8
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
        k = 2 * np.pi * nu / c                  # 书 35 行
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
    return I2, x2, (nu_arr, S_arr, dnu, c, lam0, f, w, dels)


def parabolic_extrema(p):
    """一维抛物线插值细化局部极值:返回(极大值数组, 极小值数组)。

    采样条纹直接取离散极值会低估峰、高估谷(10 样本/周期时最差 ~5%),
    用极值点处三点抛物线顶点修正到亚样本精度。
    """
    imax = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1
    imin = np.where((p[1:-1] < p[:-2]) & (p[1:-1] < p[2:]))[0] + 1
    def refine(idx):
        if len(idx) == 0:
            return np.empty(0)
        a, b, c = p[idx - 1], p[idx], p[idx + 1]
        denom = a - 2 * b + c
        delta = np.where(denom != 0, (a - c) / (2 * denom), 0.0)
        return b - (a - c) * delta / 4.0
    return refine(imax), refine(imin)


if __name__ == "__main__":
    # ----- 参数与相干量(书 9.1.1 节正文与 9.1.3 节数字) -----
    lam0, c = 650e-9, 3e8
    nu0 = c / lam0
    delnu = 2e9
    tau_c = 0.664 / delnu                       # 式 (9.6)
    l_c = c * tau_c                             # 式 (9.7)
    print(f"tau_c = 0.664/dnu = {tau_c:.3e} s(书 3.32e-10);"
          f"l_c = c*tau_c = {l_c:.4f} m(书 10 cm)")
    tau = 5e-2 / c                              # tau = deld/c = 1.67e-10 s
    g_abs = abs(gamma_temporal(tau, delnu, nu0))
    print(f"|gamma(deld/c)| 解析 = {g_abs:.4f}(书 0.6733)")

    I2, x2, (nu_arr, S_arr, dnu, c_, lam0_, f_, w_, dels_) = pc_temp()

    # ----- 数值验证 -----
    print(f"谱线权重和 sum(S*dnu) = {np.sum(S_arr) * dnu:.6f}(归一化线型应 = 1)")
    gam_num = np.sum(S_arr * np.exp(-2j * np.pi * nu_arr * tau)) * dnu
    print(f"离散度和 |sum S e^-j2pi nu tau dnu| = {abs(gam_num):.4f} "
          f"(解析 {g_abs:.4f},相对差 {abs(abs(gam_num) - g_abs) / g_abs:.2e})")

    # 解析式 (9.10) 剖面(y2=0)逐点对比
    I_an = analytic_two_beam(x2, 0.0, lam0_, f_, w_, dels_, 5e-2, c_, g_abs)
    prof = I2[I2.shape[0] // 2, :]              # 书 56 行:I2(M/2+1,:)
    rel = np.max(np.abs(prof - I_an)) / I_an.max()
    i_bad = int(np.argmax(np.abs(prof - I_an)))
    env_pk = 2 * np.pi ** 2 * w_ ** 4 / (lam0_ * f_) ** 2   # 式(9.10) 包络中心值
    print(f"剖面峰值:仿真 {prof.max():.1f},解析 {I_an.max():.1f}"
          f"(书图 9.5(b) 峰值 ~1230;包络中心值 2*pi^2*w^4/(lam0 f)^2 = {env_pk:.0f},"
          f"乘 (1+|gamma|) 的条纹峰 = {env_pk * (1 + g_abs):.0f})")
    print(f"剖面 vs 式 (9.10):最大相对偏差 = {rel:.2e}(峰值归一),"
          f"最大点在 x2 = {x2[i_bad] * 1e6:.1f} um(旁瓣区,像元化圆盘边缘采样,"
          f"即书图 9.6 评注的 lumpiness)")

    # gamma 的最小二乘恢复:prof = |G|^2 [1 + gamma cos(phi)] 对 gamma 线性
    # |G|^2 用单束圆盘的一次 FFT 直接算(与仿真同一像元化圆盘,自洽)
    per = f_ * lam0_ / dels_                    # 条纹周期 lam0 f / dels
    x1c, _ = coords(50e-3, I2.shape[0])         # 空间轴是 coords 的第一个输出!
    x1c = np.fft.fftshift(x1c)
    X1c, Y1c = np.meshgrid(x1c, x1c)
    G1 = np.fft.fft2(circle(X1c, Y1c, w_)) * (50e-3 / I2.shape[0]) ** 2 \
        / (lam0_ * f_)
    G2 = np.abs(np.fft.fftshift(G1)) ** 2       # 单束包络(居中排列,与 x2c 对齐)
    _, fx2c = coords(50e-3, I2.shape[0])
    x2c = fx2c * lam0_ * f_
    phi = 2 * np.pi * dels_ * x2c / (lam0_ * f_) + 2 * np.pi / lam0_ * 5e-2
    g2row = G2[G2.shape[0] // 2, :]             # 单束包络的 y2=0 剖面
    win = np.abs(x2c) < 2.5 * per               # 中心旁瓣区(避开第一个 Airy 零)
    A = np.stack([g2row[win], g2row[win] * np.cos(phi[win])], axis=1)
    coef, *_ = np.linalg.lstsq(A, prof[win], rcond=None)
    print(f"gamma 最小二乘恢复 = {coef[1] / coef[0]:.4f}(解析 {g_abs:.4f})")

    # 经典极值法可见度(受包络衰减 + 10 样本/周期定位偏差影响,只作量级对照)
    peak, valley = parabolic_extrema(prof[np.abs(x2) < per])
    k = min(len(peak), len(valley))
    V = np.mean((peak[:k] - valley[:k]) / (peak[:k] + valley[:k]))
    print(f"条纹可见度(中心 1 个周期插值极值)= {V:.4f}"
          f"(包络一个周期内衰减 ~5%,故比 |gamma| 略偏高)")

    # ----- 出图 1:主例子(书图 9.5) + |gamma|(图 9.2) + 采样谱(图 9.3) -----
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.6))
    ext = [x2[0], x2[-1], x2[0], x2[-1]]
    im = axes[0].imshow(I2, origin="lower", extent=ext, cmap="gray")
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
    for j, dd in enumerate([0.0, 5e-2, 10e-2, 50e-2]):
        Ij, xj, _ = pc_temp(deld=dd)
        g_j = np.exp(-(np.pi * delnu * (dd / c) / (2 * np.sqrt(np.log(2)))) ** 2)
        ext = [xj[0], xj[-1], xj[0], xj[-1]]
        axes[0, j].imshow(Ij, origin="lower", extent=ext, cmap="gray")
        axes[0, j].set_title(f"$\\Delta d$={dd*1e2:.0f} cm, "
                             f"|gamma|={g_j:.2f}")
        axes[1, j].plot(xj * 1e6, Ij[Ij.shape[0] // 2, :])
        axes[1, j].set_xlim(-200, 200)          # 书图 9.6:放大到中央 +-2e-4 m
        axes[1, j].set_title("x profile")
        axes[1, j].set_xlabel("x (um)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "out_fig9_temporal_series.png"), dpi=150)
    print("saved out_fig9_temporal.png / out_fig9_temporal_series.png")
