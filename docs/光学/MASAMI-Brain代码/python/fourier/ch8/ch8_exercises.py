# -*- coding: utf-8 -*-
"""书 8.8 节习题数值验证:8.2 离焦补偿 / 8.3 通过焦PSF / 8.5 临界采样离焦 /
8.6 Zernike vs Seidel / 8.8 Strehl 比值。

每题独立成函数,打印数值结论并保存 out_fig8_ex*.png。
依赖:python/fourier/common/base.py(coords/circle/rect)与 ch8_aberrations。
"""
import os
import sys

# 控制台打印约定:Windows 控制台默认 GBK,上标/希腊字母以外的符号
#(²、ρ 等)不在 GBK 内,日志一律用 ASCII 写法(^2、r);图形内不受限。

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, HERE)
from base import coords, circle, rect                      # noqa: E402
from ch8_aberrations import (seidel_5, psf_from_H,         # noqa: E402
                             mtf_from_psf, zernike_11)

OUT = HERE


def grids(M=1024, L=1e-3):
    """像面频率网格 + 归一化瞳坐标网格(f/5 平凸透镜参数,书 8.5 节)。"""
    lam = 0.55e-6
    wxp, zxp = 0.01, 0.1
    _, fu = coords(L, M)
    FU, FV = np.meshgrid(fu, fu)
    xhat = -lam * zxp * FU / wxp
    yhat = -lam * zxp * FV / wxp
    return dict(lam=lam, wxp=wxp, zxp=zxp, fu=fu, FU=FU, FV=FV,
                xhat=xhat, yhat=yhat, k=2 * np.pi / lam,
                f0=wxp / (lam * zxp), L=L, M=M)


def psf_of_W(g, W):
    """给波前 OPD 网格,返回 PSF(式 8.3/8.7)。"""
    H = circle(g["FU"], g["FV"], g["f0"]) * np.exp(-1j * g["k"] * W)
    return psf_from_H(H, pre_shifted=False)


# ===== 习题 8.2:离焦补偿球差 =====
def ex82():
    print("== 习题 8.2:W040=1λ 时用离焦 Wd 补偿 ==")
    g = grids(M=512)
    xh, yh = g["xhat"], g["yhat"]
    inside = np.hypot(xh, yh) <= 1.0                     # 单位圆瞳内
    w040 = 1.0 * g["lam"]
    Wd_list = np.linspace(-1.5, 0.0, 301) * g["lam"]
    rms = []
    for wd in Wd_list:
        W = seidel_5(0, 0, xh, yh, wd, w040, 0, 0, 0, 0)
        rms.append(np.sqrt(np.mean(W[inside] ** 2) - np.mean(W[inside]) ** 2))
    rms = np.array(rms)
    wd_num = Wd_list[np.argmin(rms)] / g["lam"]
    # 理论:var(c) = c^2/12 + c/6 + 4/45(c = Wd/W040)在 c = -1 取最小,
    # 即经典平衡形式 W040(r^4 - r^2)(∝ Zernike 第 11 项),σ_min = W040/sqrt(180)
    print(f"数值最优 Wd = {wd_num:+.4f} lam(理论 -W040 = -1 lam,"
          f"即 W = W040(r^4 - r^2)),最小 rms = {rms.min()/g['lam']:.4f} lam"
          f"(理论 1/sqrt(180) = {1/np.sqrt(180):.4f} lam)")
    # 无补偿的方差
    W = w040 * (np.hypot(xh, yh) ** 2) ** 2              # 纯球差
    var = np.mean(W[inside] ** 2) - np.mean(W[inside]) ** 2
    print(f"纯球差方差 = {var/g['lam']**2:.4f} lam^2(理论 4/45 = "
          f"{(1/12 + 1/180):.4f} lam^2,即书习题 8.8(d) 的 (1/12+1/180)W040^2)")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, wd, tt in [(axes[0], 0.0, "W040 only"),
                       (axes[1], -1.0 * g["lam"], "W040 + Wd=-1$\\lambda$ (balanced)")]:
        Wp = seidel_5(0, 0, xh, yh, wd, w040, 0, 0, 0, 0)
        Wp[~inside] = np.nan
        im = ax.imshow(Wp / g["lam"], origin="lower",
                       extent=[-1, 1, -1, 1], cmap="viridis")
        ax.set_title(tt + f"  (rms={np.nanstd(Wp/g['lam']):.3f}λ)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "out_fig8_ex82.png"), dpi=150)
    print("saved out_fig8_ex82.png")


# ===== 习题 8.3:通过焦 PSF(delta_d = 8(f/#)^2 Wd) =====
def ex83():
    print("== 习题 8.3:f/5 透镜轴上 PSF through-focus ==")
    g = grids(M=1024)
    fnum = g["zxp"] / (2 * g["wxp"])
    u_img, _ = coords(g["L"], g["M"])                    # 像面空间坐标(m)
    peaks, widths = [], []
    dds = [-1.0, -0.5, 0.0, 0.5, 1.0]                    # mm
    fig, axes = plt.subplots(1, 5, figsize=(15, 3.2), sharey=True)
    for ax, dd_mm in zip(axes, dds):
        Wd = dd_mm * 1e-3 / (8 * fnum ** 2)              # 逆用书式:Wd = δd/(8(f/#)^2)
        h2 = psf_of_W(g, seidel_5(0, 0, g["xhat"], g["yhat"],
                                  Wd, 4.963 * g["lam"], 0, 0, 0, 0))
        c = g["M"] // 2
        prof = h2[c, :]
        half = prof.max() / 2
        fwhm = (prof > half).sum() * (g["L"] / g["M"]) * 1e6   # um
        peaks.append(prof.max()); widths.append(fwhm)
        ax.plot(u_img * 1e3, prof / prof.max())          # 横轴:像面 u(mm)
        ax.set_xlim(-0.03, 0.03)
        ax.set_xlabel("u (mm)")
        ax.set_title(f"dd={dd_mm:+.1f} mm\nWd={Wd/g['lam']:+.2f} lam "
                     f"FWHM={fwhm:.1f} um", fontsize=9)
    print("  dd(mm)   Wd(lam)   peak        FWHM(um)")
    for dd_mm, pk, fw in zip(dds, peaks, widths):
        print(f"  {dd_mm:+.2f}   {dd_mm*1e-3/(8*fnum**2)/g['lam']:+8.2f}  "
              f"{pk:.3e}  {fw:6.1f}")
    best = int(np.argmax(peaks))
    wd_balance = -4.963                                   # 平衡球差所需 Wd(λ)
    dd_balance = 8 * fnum ** 2 * wd_balance * g["lam"] * 1e3
    print(f"最窄/最高 PSF:dd = {dds[best]:+.1f} mm(Wd = "
          f"{dds[best]*1e-3/(8*fnum**2)/g['lam']:+.2f} lam)")
    print(f"物理:离焦平衡球差要求 Wd = -W040 = -4.963 lam,即 "
          f"dd = 8(f/#)^2*Wd = {dd_balance:+.2f} mm —— -0.5 mm 档最接近,"
          f"故最锐(不是 0!这是本题的要点)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "out_fig8_ex83.png"), dpi=150)
    print("saved out_fig8_ex83.png")


# ===== 习题 8.5:最大像高点的临界采样离焦 =====
def ex85():
    print("== 习题 8.5:最大像高点 (0,1) 加正离焦到式(8.12) 临界 ==")
    g = grids(M=1024)
    L, fnum = g["L"], g["zxp"] / (2 * g["wxp"])
    lam = g["lam"]
    cf = (0.0, 4.963 * lam, 2.637 * lam, 9.025 * lam, 7.536 * lam, 0.157 * lam)
    bound = L / (4 * fnum)                               # 式 (8.12) 右端(米)
    slope0 = (4 * cf[1] + 3 * cf[2] + 2 * cf[3] + 2 * cf[4] + cf[5])  # 式 (8.16),u0=1
    wd_crit = (bound - slope0) / 2                       # 2Wd 项
    print(f"式(8.12) 界 L/(4f#) = {bound*1e6:.2f} um = {bound/lam:.2f} lam;"
          f"现有斜率(u0=1) = {slope0/lam:.2f} lam")
    print(f"临界 Wd = (界-斜率)/2 = {wd_crit/lam:+.2f} lam")
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.6))
    for ax, wdm, tt in [(axes[0], 0.5 * wd_crit, "under (Wd/2)"),
                        (axes[1], wd_crit, "critical"),
                        (axes[2], 1.5 * wd_crit, "over (1.5 Wd)")]:
        h2 = psf_of_W(g, seidel_5(0, 1, g["xhat"], g["yhat"],
                                  wdm, *cf[1:]))
        ax.imshow(np.sqrt(h2), origin="lower", cmap="gray",
                  extent=[-0.5, 0.5, -0.5, 0.5], vmax=np.sqrt(h2).max() * 0.08)
        ax.set_title(f"{tt}: Wd={wdm/lam:+.1f} lam peak={h2.max():.2e}",
                     fontsize=9)
        print(f"  Wd = {wdm/lam:+7.2f} lam -> PSF peak {h2.max():.4e}")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "out_fig8_ex85.png"), dpi=150)
    print("saved out_fig8_ex85.png")


# ===== 习题 8.6:Zernike(表 8.4 前 11 项)vs Seidel PSF/MTF =====
def ex86():
    print("== 习题 8.6:Zernike 前 11 项 vs Seidel(最大像高点)==")
    g = grids(M=1024)
    lam = g["lam"]
    # 表 8.4:ZEMAX standard Zernike 系数(单位:波长;x 0.55 μm)
    zc = np.array([7.70185740, 0.0, 1.04836665, 4.95277338, 0.0,
                   -1.83668298, 0.38477648, 0.0, -0.00829929, 0.0,
                   0.39476169]) * lam
    Wz = zernike_11(g["xhat"], g["yhat"], zc)
    cf = (0.0, 4.963 * lam, 2.637 * lam, 9.025 * lam, 7.536 * lam, 0.157 * lam)
    Ws = seidel_5(0, 1, g["xhat"], g["yhat"], *cf)
    h2z = psf_of_W(g, Wz); h2s = psf_of_W(g, Ws)
    MTFz, MTFs = mtf_from_psf(h2z), mtf_from_psf(h2s)
    print(f"PSF 峰值:Zernike {h2z.max():.4e}, Seidel {h2s.max():.4e}"
          f"(比值 {h2z.max()/h2s.max():.3f})")
    print("注:表 8.4 第 4 项含 4.95 lam 的'最佳离焦'(ZEMAX 拟合把补偿离焦"
          "并入),Seidel 组 Wd=0 —— 这是两版 PSF 不同的主要原因")
    # Zernike 拟合含 4.95λ "最佳离焦"(表 8.4 第 4 项),Seidel 的 Wd=0 ——
    # 这是两者 PSF 不同的主要原因:Zernike 版处于近"最佳焦点"位置
    fhalf = g["fu"] >= 0
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    for ax, h2, tt in [(axes[0], h2z, "Zernike PSF"),
                       (axes[1], h2s, "Seidel PSF")]:
        ax.imshow(np.sqrt(h2), origin="lower", cmap="gray",
                  extent=[-0.5, 0.5, -0.5, 0.5], vmax=np.sqrt(h2).max() * 0.08)
        ax.set_title(tt, fontsize=10)
    ax = axes[2]
    ax.plot(g["fu"][fhalf] * 1e-3, MTFz[g["M"] // 2, fhalf], label="Zernike u")
    ax.plot(g["fu"][fhalf] * 1e-3, MTFs[g["M"] // 2, fhalf], "--",
            label="Seidel u")
    ax.set_xlim(0, 40); ax.set_ylim(0, 1)
    ax.set_xlabel("f (cyc/mm)"); ax.set_ylabel("Modulation")
    ax.set_title("MTF u-axis comparison"); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "out_fig8_ex86.png"), dpi=150)
    print("saved out_fig8_ex86.png")


# ===== 习题 8.8:Strehl 比值 =====
def ex88():
    print("== 习题 8.8:Strehl 比值 vs W040 ==")
    g = grids(M=1024)
    h2_unab = psf_of_W(g, np.zeros_like(g["xhat"]))
    I0 = h2_unab.max()
    w040_waves = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 1.0])
    S_num = []
    for wv in w040_waves:
        h2 = psf_of_W(g, seidel_5(0, 0, g["xhat"], g["yhat"],
                                  wv * g["lam"], 0, 0, 0, 0, 0))
        S_num.append(h2.max() / I0)
    S_num = np.array(S_num)
    # Marechal 近似(式 8.19):S ≈ exp(-(k σ)^2),σ² = (4/45)W040²
    sig2 = (1 / 12 + 1 / 180) * (w040_waves * g["lam"]) ** 2   # 书印本漏了平方
    S_marechal = np.exp(-(g["k"] ** 2) * sig2)
    wv_fine = np.linspace(0, 1, 200)
    sig2f = (1 / 12 + 1 / 180) * (wv_fine * g["lam"]) ** 2
    print("  W040(lam)   S numeric   S Marechal")
    for wv, sn, sm in zip(w040_waves, S_num, S_marechal):
        print(f"  {wv:5.2f}      {sn:.4f}      {sm:.4f}")
    print("(e) f/5 例镜轴上(5 项 Seidel,W040=4.963 lam):见 ch8_lens_psfmtf "
          "输出 S = 0.0212;Marechal 给 exp(-86) 约 4e-38 —— 5 lam 远超 "
          "Marechal 适用域(rms < 约 0.15 lam),只能作小像差趋势参考")
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(wv_fine, np.exp(-(g["k"] ** 2) * sig2f), "-",
            label="Marechal (8.19), $\\sigma^2=\\frac{4}{45}W_{040}^2$")
    ax.plot(w040_waves, S_num, "o", label="numerical PSF peak ratio")
    ax.set_xlabel("$W_{040}$ (waves)"); ax.set_ylabel("Strehl ratio S")
    ax.set_yscale("log"); ax.set_ylim(1e-4, 1.5)
    ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "out_fig8_ex88.png"), dpi=150)
    print("saved out_fig8_ex88.png")


# ===== 习题 8.4:f/5 抛物面镜(表 8.3) =====
def ex84():
    print("== 习题 8.4:f/5 抛物面镜 PSF/MTF(lambda=0.633um,M=1024,L=0.1mm)==")
    lam = 0.633e-6
    wxp, zxp = 0.02, 0.2                       # 直径 40 mm,f = 200 mm,f/5
    L, M = 0.1e-3, 1024
    _, fu = coords(L, M)
    FU, FV = np.meshgrid(fu, fu)
    xhat = -lam * zxp * FU / wxp
    yhat = -lam * zxp * FV / wxp
    k = 2 * np.pi / lam
    f0 = wxp / (lam * zxp)
    cf = (0.0, 0.0, -1.3792 * lam, 0.4815 * lam, 0.0, 0.0)   # 表 8.3
    h2_unab = None
    for (u0, v0) in [(0, 0), (0, 1), (0.707, 0.707)]:
        W = seidel_5(u0, v0, xhat, yhat, *cf)
        H = circle(FU, FV, f0) * np.exp(-1j * k * W)
        h2 = psf_from_H(H, pre_shifted=False)
        if h2_unab is None:
            H0 = circle(FU, FV, f0)
            h2_unab = psf_from_H(H0, pre_shifted=False)
        print(f"  (u0,v0)=({u0},{v0}):  PSF peak = {h2.max():.4e}, "
              f"Strehl(对轴上衍射限) = {h2.max()/h2_unab.max():.4f}")


# ===== 习题 8.7:单项 Seidel 像差的 PSF 图(psf_map 参数) =====
def ex87():
    print("== 习题 8.7:单项像差的 PSF 图(f/20,psf_map 网格)==")
    lam = 0.5e-6
    wxp, zxp = 2.5e-3, 0.1
    fnum = zxp / (2 * wxp)
    L, M = 1e-3, 250
    k = 2 * np.pi / lam
    _, fu = coords(L, M)
    fu = np.fft.fftshift(fu)
    FU, FV = np.meshgrid(fu, fu)
    f0 = wxp / (lam * zxp)
    cases = [("1 lam defocus", (1 * lam, 0, 0, 0, 0, 0)),
             ("1 lam spherical", (0, 1 * lam, 0, 0, 0, 0)),
             ("2 lam coma", (0, 0, 2 * lam, 0, 0, 0)),
             ("1 lam astigmatism", (0, 0, 0, 1 * lam, 0, 0)),
             ("1 lam field curvature", (0, 0, 0, 0, 1 * lam, 0)),
             ("+3 lam distortion", (0, 0, 0, 0, 0, 3 * lam))]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.6))
    for ax, (tt, cf) in zip(axes.ravel(), cases):
        I = np.zeros((M, M))
        for u0 in np.arange(-0.7, 0.7001, 0.7 / 3):
            for v0 in np.arange(-0.7, 0.7001, 0.7 / 3):
                W = seidel_5(u0, v0, -2 * lam * fnum * FU,
                             -2 * lam * fnum * FV, *cf)
                H = circle(FU, FV, f0) * np.exp(-1j * k * W)
                h2 = psf_from_H(H, pre_shifted=True)
                h2 = np.roll(h2, (int(np.round(v0 * M / 2)),
                                  int(np.round(u0 * M / 2))), axis=(0, 1))
                I += h2
        ax.imshow(I ** 0.25, origin="lower", cmap="gray",
                  extent=[-0.5, 0.5, -0.5, 0.5])
        ax.set_title(tt, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "out_fig8_ex87.png"), dpi=150)
    print("saved out_fig8_ex87.png(注意:畸变项只平移 PSF、不展宽;"
          "彗差拖尾沿 v;像散随视场变椭圆/十字)")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "82"):
        ex82()
    if which in ("all", "83"):
        ex83()
    if which in ("all", "84"):
        ex84()
    if which in ("all", "85"):
        ex85()
    if which in ("all", "86"):
        ex86()
    if which in ("all", "87"):
        ex87()
    if which in ("all", "88"):
        ex88()
