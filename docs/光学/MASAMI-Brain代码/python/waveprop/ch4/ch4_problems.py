# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 4 章 =====
原书 4.3 节习题 1-5：夫琅禾费衍射与透镜焦面的五个数值/解析对比

题目（均按 4.1 节例题的网格 N=512、L=7.5 mm、Dz=20 m、lambda=1 um 复做）：
    1) 1 mm x 1 mm 方孔：解析 U = D^2 sinc(D fx) sinc(D fy)，
       fx = x2/(lambda*Dz)（式 1.67 的方孔谱，见 ch2 笔记 2.6 节）；
    2) 双缝：两个 1 mm 方孔、间隙 0.5 mm（按"边缘间隙"理解，中心距
       1.5 mm；若按中心距 0.5 mm 两孔重叠 0.5 mm，题目失义）。
       解析 = 2 D^2 sinc(D fx) sinc(D fy) cos(2 pi fx x0)，x0=0.75 mm
       ——方孔 sinc 包络 x 双缝余弦条纹；
    3) 1 mm 方孔振幅光栅 tA = 1/2[1+cos(2 pi f0 x)] rect(x/D) rect(y/D)，
       f0 = 10/D：解析 = D^2 sinc(D fy) [1/2 sinc(D fx)
       + 1/4 sinc(D(fx-f0)) + 1/4 sinc(D(fx+f0))]（零级 + 一对一级）；
    4) 1 mm 方孔相位光栅 tA = exp(i 2 pi cos(2 pi f0 x)) rect rect：
       Jacobi-Anger 展开 e^{i z cos t} = sum_n i^n J_n(z) e^{i n t} 给出
       一串衍射级，第 n 级峰值强度 = J_n^2(2 pi) (D^2/(lambda Dz))^2，
       且 sum_n J_n^2(2 pi) = 1（纯相位光栅不吞能量）。
       带宽账本（f0 = 10/D，调制深度 2 pi rad 属于深调制）：
       J_4..J_6 都不小（0.32/0.37/0.28），sum|n|<=3 只装下 31% 能量 ->
       N=512（频带 ±34133 m^-1）会混叠 69%，N=1024 仍混叠 6.2%，
       本题取 N=2048（|n|<=13，占比 100%-6e-9）才干净；
    5) 高斯束聚焦：lambda=1 um、束腰在透镜处、w=2 cm（1/e 振幅半径
       约定，U1 = exp(-r^2/w^2)）、f_l = 1 m、无限孔径。
       解析焦面场 = pi w^2/(i lambda f_l) exp(-(pi w r2/(lambda f_l))^2)
       （高斯的 FT 还是高斯），强度 1/e^2 半径 = lambda f_l/(pi w)。
       网格单独设计：L=0.4 m、N=1024（束宽 2 cm 需要大窗，焦斑 16 um
       需要细 d2 = lambda f_l/L）。

数值验证（运行自动打印）：每题数值 vs 解析的最大偏差、能量守恒、
    习题 4 的 J_n 表与带宽混叠账本、习题 5 的焦斑半径与能量。
MATLAB 注意：Bessel 函数用 besselj(1,pi*x)（jinc）与 besselj(n,2*pi)
    （习题 4）；rect/sinc/jinc 用本地函数（见配套 .m 文件末尾）。
用法：python ch4_problems.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib + scipy（jv）；base.py（sinc/rect/ft2_centered，
    同仓库 python/fourier/common）。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import rect, sinc, ft2_centered  # noqa: E402


def fraunhofer_prop(Uin, wvl, d1, Dz):
    """原书 Listing 4.1（同 ch4_fraunhofer.py）。"""
    N = Uin.shape[0]
    k = 2 * np.pi / wvl
    fX = (np.arange(N) - N // 2) / (N * d1)
    x2, y2 = np.meshgrid(wvl * Dz * fX, wvl * Dz * fX)
    Uout = np.exp(1j * k / (2 * Dz) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * Dz) * ft2_centered(Uin, d1)
    return Uout, x2, y2


def lens_against_ft(Uin, wvl, d1, f):
    """原书 Listing 4.3（同 ch4_lenses.py）。"""
    N = Uin.shape[0]
    k = 2 * np.pi / wvl
    fX = (np.arange(N) - N // 2) / (N * d1)
    x2, y2 = np.meshgrid(wvl * f * fX, wvl * f * fX)
    Uout = np.exp(1j * k / (2 * f) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * f) * ft2_centered(Uin, d1)
    return Uout, x2, y2


def err_report(tag, I_num, I_ana):
    """归一化形状误差（各自按峰值归一，扣除孔径离散化的标度差）。"""
    shape = np.abs(I_num / I_num.max() - I_ana / I_ana.max())
    scale = I_num.max() / I_ana.max() - 1
    print('    %s：标度差 %+.2f%%，归一化形状误差 max %.2e'
          % (tag, 100 * scale, shape.max()))
    return shape


def main(outdir=None):
    # ---- 4.1 节例题的公共网格 ----
    N = 512
    L = 7.5e-3
    d1 = L / N
    wvl = 1e-6
    k = 2 * np.pi / wvl
    Dz = 20.0
    n1 = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n1 * d1, n1 * d1)
    pref = np.exp(1j * k / (2 * Dz) * (x1 ** 2 * 0))  # 占位避免未用告警
    del pref

    def slice_xy(x2, U):
        """取 y2=0 行与对应 x2 升序切片。"""
        sl = U[N // 2, :]
        return x2[N // 2, :], sl

    print('===== 习题 1：1 mm 方孔 =====')
    D = 1e-3
    Uin = rect(x1 / D) * rect(y1 / D)
    Uout, x2, y2 = fraunhofer_prop(Uin, wvl, d1, Dz)
    fx = x2 / (wvl * Dz)
    A1 = D ** 2 * sinc(D * fx) * sinc(D * y2 / (wvl * Dz))
    U_th = np.exp(1j * k / (2 * Dz) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * Dz) * A1
    xs, us = slice_xy(x2, np.abs(Uout) ** 2)
    _, at = slice_xy(x2, np.abs(U_th) ** 2)
    err_report('习题 1', np.abs(Uout) ** 2, np.abs(U_th) ** 2)
    res1 = (xs, us, at)

    print('===== 习题 2：双缝（1 mm 方孔 x2，边缘间隙 0.5 mm） =====')
    x0 = 0.75e-3
    Uin = (rect((x1 - x0) / D) + rect((x1 + x0) / D)) * rect(y1 / D)
    Uout, x2, y2 = fraunhofer_prop(Uin, wvl, d1, Dz)
    A2 = 2 * D ** 2 * sinc(D * x2 / (wvl * Dz)) \
        * sinc(D * y2 / (wvl * Dz)) * np.cos(2 * np.pi * x2 / (wvl * Dz) * x0)
    U_th = np.exp(1j * k / (2 * Dz) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * Dz) * A2
    xs, us = slice_xy(x2, np.abs(Uout) ** 2)
    _, at = slice_xy(x2, np.abs(U_th) ** 2)
    err_report('习题 2', np.abs(Uout) ** 2, np.abs(U_th) ** 2)
    res2 = (xs, us, at)

    print('===== 习题 3：振幅光栅（f0 = 10/D） =====')
    f0 = 10 / D
    Uin = 0.5 * (1 + np.cos(2 * np.pi * f0 * x1)) * rect(x1 / D) * rect(y1 / D)
    Uout, x2, y2 = fraunhofer_prop(Uin, wvl, d1, Dz)
    fxg = x2 / (wvl * Dz)
    A3 = D ** 2 * sinc(D * y2 / (wvl * Dz)) * (
        0.5 * sinc(D * fxg)
        + 0.25 * sinc(D * (fxg - f0))
        + 0.25 * sinc(D * (fxg + f0)))
    U_th = np.exp(1j * k / (2 * Dz) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * Dz) * A3
    xs, us = slice_xy(x2, np.abs(Uout) ** 2)
    _, at = slice_xy(x2, np.abs(U_th) ** 2)
    err_report('习题 3', np.abs(Uout) ** 2, np.abs(U_th) ** 2)
    res3 = (xs, us, at)

    print('===== 习题 4：相位光栅（f0 = 10/D，N=2048） =====')
    N4 = 2048
    d14 = L / N4
    n4 = np.arange(N4) - N4 // 2
    x1b, y1b = np.meshgrid(n4 * d14, n4 * d14)
    Uin = np.exp(2j * np.pi * np.cos(2 * np.pi * f0 * x1b)) \
        * rect(x1b / D) * rect(y1b / D)
    Uout, x2b, y2b = fraunhofer_prop(Uin, wvl, d14, Dz)
    I4 = np.abs(Uout) ** 2
    fx4 = x2b / (wvl * Dz)
    nn = np.arange(-8, 9)
    A4 = np.zeros_like(x2b, dtype=complex)
    for n_ord in nn:
        A4 += 1j ** n_ord * jv(n_ord, 2 * np.pi) * D ** 2 \
            * sinc(D * (x2b / (wvl * Dz) - n_ord * f0)) \
            * sinc(D * y2b / (wvl * Dz))
    U_th4 = np.exp(1j * k / (2 * Dz) * (x2b ** 2 + y2b ** 2)) \
        / (1j * wvl * Dz) * A4
    err_report('习题 4', I4, np.abs(U_th4) ** 2)
    scale_pk = (D ** 2 / (wvl * Dz)) ** 2
    print('    各衍射级峰值辐照度（理论 = J_n^2(2pi) x %.4e）：' % scale_pk)
    row = I4[N4 // 2, :]
    fx_row = x2b[N4 // 2, :] / (wvl * Dz)
    for n_ord in [0, 1, 2, 3]:
        j = np.argmin(np.abs(fx_row - n_ord * f0))
        print('      n=%+d：数值 %.4e vs 理论 %.4e（J_%d(2pi)=%.4f）'
              % (n_ord, row[j], jv(n_ord, 2 * np.pi) ** 2 * scale_pk,
                 n_ord, jv(n_ord, 2 * np.pi)))
    s = sum(jv(n_ord, 2 * np.pi) ** 2 for n_ord in range(-40, 41))
    print('    sum_n J_n^2(2pi) = %.6f（=1，纯相位光栅不吞能量）' % s)
    s_all = sum(jv(n_ord, 2 * np.pi) ** 2 for n_ord in range(-60, 61))
    for tag, Nq in [('N=512（原书网格）', 512), ('N=1024', 1024),
                    ('N=2048（本题）', 2048)]:
        nb = int(1 / (2 * (L / Nq) * f0))       # 频带内能装下的最高级数
        s_in = sum(jv(n_ord, 2 * np.pi) ** 2 for n_ord in range(-nb, nb + 1))
        print('    %s：频带 ±%.0f m^-1 装到 |n|<=%d，能量占比 %.1f%%'
              '（混叠 %.2f%%）'
              % (tag, 1 / (2 * (L / Nq)), nb,
                 100 * s_in / s_all, 100 * (1 - s_in / s_all)))
    P4 = np.sum(I4) * (wvl * Dz / L) ** 2
    print('    能量：sum|U2|^2 d2^2 = %.4e vs 物面积 D^2 = %.4e（相对偏差 %.2e）'
          % (P4, D ** 2, abs(P4 - D ** 2) / D ** 2))
    res4 = (x2b[N4 // 2, :], row,
            np.abs(U_th4[N4 // 2, :]) ** 2)

    print('===== 习题 5：高斯束聚焦（w=2 cm，f_l=1 m，无限孔径） =====')
    w = 2e-2
    fl = 1.0
    L5, N5 = 0.4, 1024
    d15 = L5 / N5
    n5 = np.arange(N5) - N5 // 2
    x1c, y1c = np.meshgrid(n5 * d15, n5 * d15)
    Uin5 = np.exp(-(x1c ** 2 + y1c ** 2) / w ** 2)
    Uout5, x2c, y2c = lens_against_ft(Uin5, wvl, d15, fl)
    r2 = np.hypot(x2c, y2c)
    U_th5 = np.pi * w ** 2 / (1j * wvl * fl) \
        * np.exp(-(np.pi * w * r2 / (wvl * fl)) ** 2)
    I5n, I5a = np.abs(Uout5) ** 2, np.abs(U_th5) ** 2
    shape5 = np.abs(I5n / I5n.max() - I5a / I5a.max())
    print('    采样：d1=%.3e m（每 w 跨 %.0f 个样本），d2=%.1e m'
          '（焦斑 1/e^2 半径 %.2e m 跨 %.1f 个样本）'
          % (d15, w / d15, wvl * fl / L5, wvl * fl / (np.pi * w),
             wvl * fl / (np.pi * w) / (wvl * fl / L5)))
    print('    峰值辐照度：数值 %.4e vs 理论 (pi w^2/(lambda f_l))^2 = %.4e'
          % (I5n.max(), (np.pi * w ** 2 / (wvl * fl)) ** 2))
    print('    归一化形状误差 max = %.2e' % shape5.max())
    row5 = I5n[N5 // 2, :]
    x5row = x2c[N5 // 2, :]
    half = row5.max() / np.e ** 2
    i_pk = row5.argmax()
    right = i_pk + np.argmax(row5[i_pk:] < half)
    # 阈值穿点线性内插（相邻样本 2.5 um，直接取样本会有 +0.6 样本量化差）
    xa, xb2 = x5row[right - 1], x5row[right]
    ya, yb2 = row5[right - 1], row5[right]
    x_cross = xb2 - (yb2 - half) / (yb2 - ya) * (xb2 - xa)
    print('    强度 1/e^2 半径：数值 %.4e m vs 理论 lambda f_l/(pi w) = %.4e m'
          % (x_cross - x5row[i_pk], wvl * fl / (np.pi * w)))
    P5 = np.sum(I5n) * (wvl * fl / L5) ** 2
    P5th = np.pi * w ** 2 / 2
    print('    能量：sum|U2|^2 d2^2 = %.4e vs 束总功率 pi w^2/2 = %.4e'
          '（相对偏差 %.2e）' % (P5, P5th, abs(P5 - P5th) / P5th))
    res5 = (x2c[N5 // 2, :], row5, I5a[N5 // 2, :])

    if outdir:
        # ---- 图 fig-4-8：习题 1-4（习题 3/4 的 ±1 级在 x2=±lambda*Dz*f0
        #      = ±200 mm，显示范围放宽才能看到级次结构） ----
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6))
        titles = ['习题 1：1 mm 方孔',
                  '习题 2：双缝（间隙 0.5 mm）',
                  '习题 3：振幅光栅 $f_0=10/D$',
                  '习题 4：相位光栅 $f_0=10/D$（N=2048）']
        halves = [0.1, 0.1, 0.3, 0.7]
        for ax, res, ti, hw in zip(axes.ravel(),
                                   [res1, res2, res3, res4], titles, halves):
            xs, us, at = res
            m = np.abs(xs) <= hw
            ax.plot(xs[m] * 1e3, at[m], color='0.6', lw=1.6, ls='--',
                    label='解析')
            ax.plot(xs[m] * 1e3, us[m], 'b.', ms=3, label='数值')
            ax.set_xlabel(r'$x_2$ [mm]')
            ax.set_ylabel(r'辐照度 [W/m$^2$]')
            ax.set_title(ti, fontsize=11)
            ax.legend(fontsize=9)
            ax.grid(alpha=0.25)
        fig.tight_layout()
        fp = os.path.join(outdir, 'fig-4-8.png')
        fig.savefig(fp, dpi=150)
        print('\n已保存：%s' % fp)

        # ---- 图 fig-4-9：习题 5 ----
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
        msrc = np.abs(x1c[0, :]) <= 0.05
        axes[0].plot(x1c[0, :][msrc], Uin5[N5 // 2, :][msrc], 'k-', lw=1.4)
        axes[0].axvline(-w, color='r', ls=':', lw=1)
        axes[0].axvline(w, color='r', ls=':', lw=1)
        axes[0].set_xlabel(r'$x_1$ [m]')
        axes[0].set_ylabel('振幅')
        axes[0].set_title('(a) 透镜面振幅（红线 $\\pm w$，窗 L=0.4 m 截断'
                          ' $e^{-100}$）')
        mf = np.abs(res5[0]) <= 6e-5
        axes[1].plot(res5[0][mf] * 1e6, res5[2][mf], color='0.6', lw=1.8,
                     ls='--', label='解析')
        axes[1].plot(res5[0][mf] * 1e6, res5[1][mf], 'b.', ms=4, label='数值')
        axes[1].set_xlabel(r'$x_2$ [$\mu$m]')
        axes[1].set_ylabel(r'辐照度 [W/m$^2$]')
        axes[1].set_title('(b) 焦面光斑（1/e$^2$ 半径 15.9 $\\mu$m）')
        axes[1].legend(fontsize=9)
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.tight_layout()
        fp = os.path.join(outdir, 'fig-4-9.png')
        fig.savefig(fp, dpi=150)
        print('已保存：%s' % fp)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
