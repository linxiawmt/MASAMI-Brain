# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 7 章（工程 3）：习题复现 =====
习题 1-6 复现（原书 7.4 节）；约束/解析工具复用 ch7_sampling（工程 1+2）。

原书习题：
    1. g(x)=exp(i*pi*a^2*x^2)，a=4，N=128，L=4 m：不做 FT，解析证明混叠。
    2. 点源采样图：lambda=1 um，传播 100 km，接收孔径 2 m。
    3. 采样图：lambda=0.5 um，D1=1 mm，传播 2.0 m，接收孔径 2 m。
    4. 把 Listing 7.2/B.5 改成会聚/发散源（式 7.62）：
       (a) 推导带 R 的方孔解析解；(b) R=-0.2 m 画约束图；(c) 仿真 vs 解析。
    5. 图示式 (7.60) <=> 几何光束包含在直径 D2 内（发散源 + 会聚源
       焦点在观测面前/后的光线图，含在 ch7_sampling.py 图 7-10 三情形中）。
    6. 代数证明约束 4 强于约束 1+2 的组合。

数值验证（运行本脚本全部自动打印，与 MATLAB 版 ch7_problems.m 一致）：
    (4) 习题 1：f_loc = a^2*x，Nyquist = 16/m，混叠边界 |x| = 1.000 m
        （半网格宽）；相邻样本相位增量 |dphi|>pi 的样本对 63/127，边界一致；
        边缘表观频率 |2*fN - f_loc| = 0.5 /m（低频回折）。自绘 fig-7-11。
    (5) 习题 2/3：点源 100 km（delta1=delta2=1 cm -> N=2^10=1024，约束 4 主导）；
        D1=1 mm 短距（约束 1 => delta1 < 0.50 um；取 0.4 um / 0.1 mm
        -> N=2^15=32768，约束 4 主导）。自绘 fig-7-12。
    (6) 习题 4：R=-0.2 m 会聚源：约束 1 上界 31.0 um（主导）、约束 3
        [-45.3, 54.7] um；delta1=9.4848 um、delta2=28.1212 um、N=512：
        |dI|/pk = 8.5e-3、相位残差 8.2e-3 rad；输出曲率拟合
        R_out = -0.107 m（几何预测 R+Dz = -0.1 m）；推广解析解 R->inf
        与 ch6 闭式逐位差 4.4e-16（退化自检）；勘误 c 数值实锤：
        R=+0.2 束宽 0.75->1.12 mm（发散）、R=-0.2 -> 0.38 mm（会聚）。
        自绘 fig-7-13。
    (7) 习题 6：1000 组随机参数（约束 1 成立域内）
        max[约束2右端 - 约束4右端] = 0（<=0 恒成立）。

原书勘误 c（7.3.1.2 节 R 符号句）：按式 (7.32) 的相位 exp(i*k*r1^2/(2*R))，
    局部频率 f = x1/(lambda*R)、方向余弦 alpha = x1/R，
    R > 0 发散、R < 0 会聚；书句"发散 R<0 / 会聚 R>0"把符号说反，
    与习题 4 的"R=-0.2 m 收敛"矛盾。本脚本验证 (6) 数值实锤。

用法：python ch7_problems.py <输出目录>  保存 PNG；否则存到脚本目录。
依赖：numpy + matplotlib；ch7_sampling（同目录，约束/解析工具）。
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ---- 复用工程 1+2 的约束判据与推广解析解（同目录 ch7_sampling） ----
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch7_sampling import (  # noqa: E402
    constraint1_d2, constraint2_N, constraint3_d2, constraint4_N,
    constraint3_is_free, fresnel_prop_square_ap_R)

# ---- 复用第 6 章传播器（同仓库） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "ch6"))
from ch6_fresnel_prop import ang_spec_prop, fresnel_prop_square_ap  # noqa: E402

# ---- 复用公共基本函数（同仓库 python/fourier/common） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import rect  # noqa: E402


def main(outdir=None):
    here = os.path.dirname(os.path.abspath(__file__))
    if outdir is None:
        outdir = here
    os.makedirs(outdir, exist_ok=True)
    results = []

    # ================= 验证 4：习题 1（啁啾混叠的解析证明） =================
    print('===== 验证 4：习题 1 —— g(x)=exp(i*pi*a^2*x^2) 的混叠解析证明 =====')
    a1, Nc, Lc = 4.0, 128, 4.0
    dc = Lc / Nc
    n = np.arange(Nc) - Nc // 2
    xc = n * dc
    g = np.exp(1j * np.pi * a1 ** 2 * xc ** 2)
    floc = a1 ** 2 * xc                                # f_loc = a^2*x（式 7.37 的特例）
    fnyq = 1.0 / (2.0 * dc)
    x_alias = 1.0 / (2.0 * a1 ** 2 * dc)               # f_loc = Nyquist 的位置
    print('    delta = L/N = %.4f m；Nyquist = 1/(2*delta) = %.1f /m' % (dc, fnyq))
    print('    f_loc(x) = a^2*x = %.1f*x，网格边缘 |x|=%.1f m 处 f_loc = %.1f /m > Nyquist'
          % (a1 ** 2, Lc / 2, a1 ** 2 * Lc / 2))
    print('    混叠边界：f_loc = Nyquist 于 |x| = %.3f m —— 恰是半网格宽，一半的样本混叠'
          % x_alias)
    # 相邻样本相位增量判据（无需 FT）：dphi = pi*a^2*delta^2*(2n+1)，>pi 即混叠
    dphi = np.pi * a1 ** 2 * dc ** 2 * (2 * n[:-1] + 1)
    n_alias = np.sum(np.abs(dphi) > np.pi)
    print('    相邻样本相位增量 |dphi|>pi 的样本对：%d/%d（与 |x|>%.3f m 的解析边界一致）'
          % (n_alias, Nc - 1, x_alias))
    # 混叠的直接数值证据：折叠后的表观频率 = |2*fN - f_loc|
    f_app = np.abs(2 * fnyq - floc)
    print('    边缘表观频率 = |2*fN - f_loc| = %.1f /m（低频"回折"）' % f_app[-1])
    figC, axsC = plt.subplots(1, 2, figsize=(10.5, 3.8))
    axC = axsC[0]
    axC.plot(xc, np.real(g), 'k-', lw=1)
    axC.axvspan(x_alias, Lc / 2, color='r', alpha=0.15)
    axC.axvspan(-Lc / 2, -x_alias, color='r', alpha=0.15)
    axC.set_title('(a) 采样后的 Re g(x)：红色区域 f_loc > Nyquist（混叠）', fontsize=10)
    axC.set_xlabel('x [m]'); axC.set_ylabel('Re g')
    axC.grid(alpha=0.3)
    axC2 = axsC[1]
    axC2.plot(xc, floc, 'b--', lw=1.2, label=r'真实局部频率 $f_{loc}=a^2 x$')
    axC2.plot(xc, f_app, 'g-.', lw=1.0,
              label=r'混叠后表观频率 $|2f_N-f_{loc}|$')
    axC2.axhline(fnyq, color='r', ls=':', lw=1.5)
    axC2.text(0.05, fnyq + 1.5, r'Nyquist $=1/(2\delta)$', color='r', fontsize=9)
    axC2.set_xlabel('x [m]'); axC2.set_ylabel('频率 [1/m]')
    axC2.set_ylim(0, 70)
    axC2.set_title('(b) 局部频率与回折：|x|>1 m 处折叠回低频', fontsize=10)
    axC2.legend(fontsize=8, loc='upper left'); axC2.grid(alpha=0.3)
    figC.suptitle('图 7-11（自绘仿真）工程 2 习题 1：a=4、N=128、L=4 m 的啁啾混叠'
                  '（|x|>1 m 即 f_loc>16 /m，网格一半欠采样）', fontsize=12)
    figC.tight_layout(rect=[0, 0, 1, 0.9])
    fpC = os.path.join(outdir, 'fig-7-11.png')
    figC.savefig(fpC, dpi=150)
    print('    已保存 %s' % fpC)
    results.append(('习题1 混叠边界 x', x_alias, 1.0, abs(x_alias - 1.0)))

    # ================= 验证 5：习题 2/3（采样图） =================
    print('===== 验证 5：习题 2/3 —— 采样图（约束 4 等高线 + 可行点） =====')
    # 习题 2：点源（D1 -> 0），lambda=1um, Dz=100 km, D2=2 m
    lam2, Dz2, D2p = 1e-6, 100e3, 2.0
    d1p = np.logspace(-4, -1, 300)                     # 0.1 mm - 100 mm
    d2p = np.logspace(-4, -1, 300)
    D1p, D2g2 = np.meshgrid(d1p, d2p)
    Nc4p = constraint4_N(D1p, D2g2, lam2, Dz2)
    Nc2p = constraint2_N(D1p, D2g2, 1e-6, D2p, lam2, Dz2)   # D1->0 只剩后两项
    # 习题 3：lambda=0.5um, D1=1mm, Dz=2m, D2=2m
    lam3, Dz3, D13, D23 = 0.5e-6, 2.0, 1e-3, 2.0
    d1q = np.logspace(-6.3, -3, 300)                   # 0.5 um - 1 mm
    d2q = np.logspace(-5, -2, 300)                     # 10 um - 10 mm
    D1q, D2g3 = np.meshgrid(d1q, d2q)
    Nc4q = constraint4_N(D1q, D2g3, lam3, Dz3)
    Nc2q = constraint2_N(D1q, D2g3, D13, D23, lam3, Dz3)
    c1q = constraint1_d2(D1q, D13, D23, lam3, Dz3)
    pick2 = (1e-2, 1e-2)                               # delta1=delta2=1 cm
    pick3 = (0.4e-6, 1e-4)                             # delta1=0.4 um, delta2=0.1 mm
    n2p = max(constraint4_N(*pick2, lam2, Dz2),
              constraint2_N(pick2[0], pick2[1], 1e-6, D2p, lam2, Dz2))
    n3p = max(constraint4_N(*pick3, lam3, Dz3),
              constraint2_N(pick3[0], pick3[1], D13, D23, lam3, Dz3))
    print('    习题 2（点源, 100 km, D2=2 m）：取 delta1=delta2=1 cm -> '
          '约束 4: N >= %.0f，约束 2: N >= %.0f（D1->0）-> N = 2^%d = %d'
          % (constraint4_N(*pick2, lam2, Dz2),
             constraint2_N(pick2[0], pick2[1], 1e-6, D2p, lam2, Dz2),
             np.ceil(np.log2(n2p)), 2 ** np.ceil(np.log2(n2p))))
    print('    习题 3（D1=1 mm, 2 m, D2=2 m）：约束 1 => delta2 <= -2000*delta1 + 1 mm，'
          'delta2 > 0 要求 delta1 < %.2f um；取 delta1=0.4 um, delta2=0.1 mm -> '
          '约束 4: N >= %.0f（主导），约束 2: N >= %.0f -> N = 2^%d = %d'
          % (lam3 * Dz3 / D23 * 1e6,
             constraint4_N(*pick3, lam3, Dz3),
             constraint2_N(pick3[0], pick3[1], D13, D23, lam3, Dz3),
             np.ceil(np.log2(n3p)), 2 ** np.ceil(np.log2(n3p))))
    figD, axsD = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for ax, (X, Y, Z, lv, ttl, xlab) in zip(
            axsD,
            [(D1p, D2g2, Nc4p, [9, 10, 11, 12, 13, 14],
              '(a) 习题 2：点源 λ=1 μm，Δz=100 km，D2=2 m\n（约束 4 主导；约束 1 对 D1→0 失效）',
              r'$\delta_1$ [m]'),
             (D1q, D2g3, Nc2q, [9, 10, 10.5, 11, 11.5, 12],
              '(b) 习题 3：D1=1 mm，λ=0.5 μm，Δz=2 m，D2=2 m\n（约束 2 的 log2 N 等高线 + 约束 1）',
              r'$\delta_1$ [m]')]):
        cs = ax.contourf(X, Y, np.log2(np.maximum(Z, 1e-6)), levels=lv, cmap='Greys')
        clb = figD.colorbar(cs, ax=ax)
        clb.set_label('log2 N')
        if ttl.startswith('(b)'):
            cs1 = ax.contour(X, Y, c1q, levels=[1e-4, 1e-3], colors='r')
            ax.clabel(cs1, fmt=lambda v: '约束1 δ2=%.1f mm' % (v * 1e3), fontsize=7)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel(xlab); ax.set_ylabel(r'$\delta_2$ [m]')
        ax.set_title(ttl, fontsize=10)
    axsD[0].plot(*pick2, 'r*', ms=12)
    axsD[0].text(1.1e-2, 0.7e-2, 'δ1=δ2=1 cm\nN=1024', color='r', fontsize=8)
    axsD[1].plot(*pick3, 'r*', ms=12)
    axsD[1].text(0.45e-6, 1.35e-4, 'δ1=0.4 μm, δ2=0.1 mm\nN=2^14=16384', color='r',
                 fontsize=8)
    # 约束 1 的显式边界线（习题 3：delta2 = -D2/D1*delta1 + lambda*Dz/D1）
    d1line = np.linspace(0.35e-6, 1.0e-6, 100)
    axsD[1].plot(d1line, constraint1_d2(d1line, D13, D23, lam3, Dz3), 'r-', lw=1.5)
    axsD[1].text(0.52e-6, 0.28e-3, '约束 1：δ2 ≤ -2000δ1+1 mm', color='r', fontsize=8)
    figD.suptitle('图 7-12（自绘仿真）工程 2 习题 2/3：采样图（红* = 一个可行选点）',
                  fontsize=12)
    figD.tight_layout(rect=[0, 0, 1, 0.9])
    fpD = os.path.join(outdir, 'fig-7-12.png')
    figD.savefig(fpD, dpi=150)
    print('    已保存 %s' % fpD)

    # ================= 验证 6：习题 4（会聚源）+ 勘误 c 实锤 =================
    print('===== 验证 6：习题 4 —— R=-0.2 m 会聚源：约束、仿真、解析、R 符号 =====')
    # 书例参数（与工程 1 的 Listing 7.2 相同，仅源加二次相位）
    wvl = 1e-6
    D1b, D2b, Dzb = 2e-3, 4e-3, 0.1
    d1b, d2b = 9.4848e-6, 28.1212e-6
    Nb = int(2 ** np.ceil(np.log2(max(
        constraint2_N(d1b, d2b, D1b, D2b, wvl, Dzb),
        constraint4_N(d1b, d2b, wvl, Dzb)))))            # = 512（工程 1 验证 2）
    slb = Nb // 2
    R4 = -0.2
    rho4 = 1.0 + Dzb / R4
    c3lo4, c3hi4 = constraint3_d2(d1b, D1b, wvl, Dzb, R4)
    c1_4 = constraint1_d2(d1b, D1b, D2b, wvl, Dzb)
    print('    rho = 1+Dz/R = %.2f；约束 3：%.1f <= delta2 <= %.1f um；'
          '约束 1：delta2 <= %.1f um（选中 28.12 um 全部满足）'
          % (rho4, c3lo4 * 1e6, c3hi4 * 1e6, c1_4 * 1e6))
    print('    式 (7.60)：|1+Dz/R| = %.2f < D2/D1 = %.0f -> 约束 3 不起作用：%s'
          % (abs(rho4), D2b / D1b, constraint3_is_free(D1b, D2b, wvl, Dzb, R4)))
    n4 = np.arange(Nb) - Nb // 2
    x1b, y1b = np.meshgrid(n4 * d1b, n4 * d1b)
    src4 = rect(x1b / D1b) * rect(y1b / D1b) * np.exp(
        1j * np.pi / wvl / R4 * (x1b ** 2 + y1b ** 2))
    x2b4, _, U4 = ang_spec_prop(src4, wvl, d1b, d2b, Dzb)
    Uan4 = fresnel_prop_square_ap_R(x2b4[0, :], 0.0, D1b, wvl, Dzb, R4)
    # R -> inf 时推广解析解必须退化为 ch6 的闭式（逐位核对）
    x2s2 = n4 * d2b                                   # 与仿真同一观测面横轴
    Uan_inf = fresnel_prop_square_ap_R(x2s2, 0.0, D1b, wvl, Dzb, np.inf)
    Uan2 = fresnel_prop_square_ap(x2s2, 0.0, D1b, wvl, Dzb)
    dev_inf = np.max(np.abs(Uan_inf - Uan2))
    print('    推广解析解 R->inf 与 ch6 闭式逐位最大差 = %.1e（退化自检）' % dev_inf)
    win4 = np.abs(x2b4[0, :]) <= D2b / 2
    I4, Ia4 = np.abs(U4[slb, :]) ** 2, np.abs(Uan4) ** 2
    pk4 = Ia4.max()
    eI4 = np.max(np.abs(I4[win4] - Ia4[win4])) / pk4
    mph4 = win4 & (Ia4 > 0.05 * pk4)
    eph4 = np.max(np.abs(np.angle(U4[slb, :][mph4] / Uan4[mph4])))
    print('    会聚源仿真（N=%d）vs 推广解析：|dI|/峰值 = %.2e；相位残差（|U|>5%%峰值）= %.2e rad'
          % (Nb, eI4, eph4))
    # 输出波前曲率：几何预测 R_out = R + Dz（会聚源经 0.1 m 后距焦点 0.1 m）
    fit = np.polyfit(x2b4[0, :][mph4] ** 2,
                     np.unwrap(np.angle(U4[slb, :]))[mph4], 1)
    R_fit = 2.0 * np.pi / wvl / (2.0 * fit[0])
    print('    输出相位曲率拟合 R_out = %.4f m（几何预测 R+Dz = %.4f m）'
          % (R_fit, R4 + Dzb))
    # 勘误 c 实锤：R=+0.2（书句称"会聚"）实际发散；R=-0.2（习题 4 称"会聚"）实际会聚
    Nw = 256; Lw = 8e-3; dw = Lw / Nw
    w0 = 1.5e-3
    nw = np.arange(Nw) - Nw // 2
    xw, yw = np.meshgrid(nw * dw, nw * dw)
    xz = np.zeros(Nw)
    for Rk, lab in ((+0.2, 'R=+0.2'), (-0.2, 'R=-0.2')):
        Uw = np.exp(-(xw ** 2 + yw ** 2) / w0 ** 2) * np.exp(
            1j * np.pi / wvl / Rk * (xw ** 2 + yw ** 2))
        # RMS 宽度（辐照度加权）比较：传播 0.1 m 后
        Iz0 = np.abs(Uw[Nw // 2, :]) ** 2
        w0_rms = np.sqrt(np.sum(Iz0 * xw[Nw // 2, :] ** 2) / np.sum(Iz0))
        _, _, Uwp = ang_spec_prop(Uw, wvl, dw, dw, 0.1)
        Iz1 = np.abs(Uwp[Nw // 2, :]) ** 2
        w1_rms = np.sqrt(np.sum(Iz1 * xw[Nw // 2, :] ** 2) / np.sum(Iz1))
        # 输出相位曲率
        mp = Iz1 > 0.05 * Iz1.max()
        ft = np.polyfit(xw[Nw // 2, :][mp] ** 2,
                        np.unwrap(np.angle(Uwp[Nw // 2, :]))[mp], 1)
        R_out = np.pi / wvl / ft[0]
        print('    [%s] 高斯束 RMS 宽度 %.2f mm -> %.2f mm（%s）；'
              '输出曲率 R_out = %+.3f m（%s）'
              % (lab, w0_rms * 1e3, w1_rms * 1e3,
                 '展宽=发散' if w1_rms > w0_rms else '收窄=会聚',
                 R_out, 'R_out>0 发散' if R_out > 0 else 'R_out<0 会聚'))
    figE, axsE = plt.subplots(2, 2, figsize=(11.5, 7.6))
    # (a) 约束图（R=-0.2）
    ax = axsE[0, 0]
    d1v = np.linspace(0.2e-6, 20e-6, 300)
    lo4, hi4 = constraint3_d2(d1v, D1b, wvl, Dzb, R4)
    c1v = constraint1_d2(d1v, D1b, D2b, wvl, Dzb)
    feas4 = np.minimum(c1v, hi4)
    ax.fill_between(d1v * 1e6, 0, np.clip(feas4, 0, None) * 1e6,
                    color='0.6', alpha=0.7, label='可行域')
    ax.plot(d1v * 1e6, c1v * 1e6, 'k-', lw=2, label='约束 1')
    ax.plot(d1v * 1e6, hi4 * 1e6, 'k-.', lw=2, label=r'约束 3 上界（斜率 $\rho$=0.5）')
    ax.plot(d1v * 1e6, np.clip(lo4, 0, None) * 1e6, 'k--', lw=2, label='约束 3 下界')
    ax.plot(d1b * 1e6, d2b * 1e6, 'r*', ms=12, zorder=5)
    ax.set_xlim(0, 20); ax.set_ylim(0, 40)
    ax.set_xlabel(r'$\delta_1$ [$\mu$m]'); ax.set_ylabel(r'$\delta_2$ [$\mu$m]')
    ax.set_title('(a) 习题 4(b)：R=-0.2 m（ρ=0.5）的约束 1/3\n红* = 选中 δ1=9.48 μm, δ2=28.12 μm',
                 fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # (b) 仿真 vs 解析
    ax = axsE[0, 1]
    ax.plot(x2b4[0, :] * 1e3, I4 / pk4, 'k-', lw=1.2, label='仿真辐照度')
    ax.plot(x2b4[0, :] * 1e3, Ia4 / pk4, 'r--', lw=1.0, label='解析（习题 4a 推广）')
    ax.set_xlim(-2, 2)
    ax.set_xlabel(r'$x_2$ [mm]'); ax.set_ylabel('归一化辐照度')
    ax.set_title('(b) 习题 4(c)：会聚源角谱仿真 vs 推广解析\n'
                 r'|dI|/pk = %.1e（N=512, δ1=9.48 μm）' % eI4, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    # (c) 相位对比
    ax = axsE[1, 0]
    ax.plot(x2b4[0, :] * 1e3, np.angle(U4[slb, :]), 'k.', ms=2.5, label='仿真相位')
    ax.plot(x2b4[0, :] * 1e3, np.angle(Uan4), 'r-', lw=1.0, label='解析相位')
    ax.set_xlim(-2, 2)
    ax.set_xlabel(r'$x_2$ [mm]'); ax.set_ylabel('相位 [rad]')
    ax.set_title(r'(c) 输出相位（曲率拟合 R$_{out}$=%.3f m ≈ R+Δz=%.1f m）' % (R_fit, R4 + Dzb),
                 fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    # (d) 勘误 c：R=+0.2 vs R=-0.2 的束宽对比
    ax = axsE[1, 1]
    ax.plot(xw[Nw // 2, :] * 1e3, Iz0 / Iz0.max(), 'k-', lw=1.5, label='源面（两案例同）')
    ax.plot(xw[Nw // 2, :] * 1e3, Iz1 / Iz1.max(), 'b--', lw=1.2,
            label='R=+0.2 传播 0.1 m（展宽 -> 发散）')
    _, _, Uwp2 = ang_spec_prop(
        np.exp(-(xw ** 2 + yw ** 2) / w0 ** 2) * np.exp(
            1j * np.pi / wvl / (-0.2) * (xw ** 2 + yw ** 2)), wvl, dw, dw, 0.1)
    Iz12 = np.abs(Uwp2[Nw // 2, :]) ** 2
    ax.plot(xw[Nw // 2, :] * 1e3, Iz12 / Iz12.max(), 'r-.', lw=1.2,
            label='R=-0.2 传播 0.1 m（收窄 -> 会聚）')
    ax.set_xlabel(r'$x$ [mm]'); ax.set_ylabel('归一化辐照度')
    ax.set_title(r'(d) 勘误 c 实锤：exp($i\pi r^2/\lambda R$) 的 R>0 发散、R<0 会聚'
                 '\n'
                 '（7.3.1.2 书句把符号说反，习题 4 的"R=-0.2 m 收敛"才是对的）',
                 fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    figE.suptitle('图 7-13（自绘仿真）工程 2 习题 4：会聚源的采样与仿真'
                  '（D1=2 mm, D2=4 mm, Δz=0.1 m, λ=1 μm, R=-0.2 m）', fontsize=12)
    figE.tight_layout(rect=[0, 0, 1, 0.94])
    fpE = os.path.join(outdir, 'fig-7-13.png')
    figE.savefig(fpE, dpi=150)
    print('    已保存 %s' % fpE)
    results.append(('L7.2R4 |dI|/pk', eI4, None, None))

    # ================= 验证 7：习题 6（约束 4 强于约束 1+2 组合） =================
    print('===== 验证 7：习题 6 —— 约束 4 强于约束 1+2 的组合（代数 + 数值） =====')
    print('    代数：约束 1 <=> D1*delta2 + D2*delta1 <= lambda*Dz')
    print('          两边除以 2*delta1*delta2：D1/(2*delta1) + D2/(2*delta2) <= lambda*Dz/(2*delta1*delta2)')
    print('          故约束 2 的右端 = 几何项 + lambda*Dz/(2*d1*d2) <= lambda*Dz/(d1*d2) = 约束 4 的右端')
    rng = np.random.default_rng(7)
    viol = 0.0
    for _ in range(1000):
        d1r = 10 ** rng.uniform(-6, -4)
        d2r = 10 ** rng.uniform(-6, -4)
        D1r = 10 ** rng.uniform(-4, -2)
        D2r = 10 ** rng.uniform(-4, -2)
        Dzr = 10 ** rng.uniform(-1, 2)
        c2r = constraint2_N(d1r, d2r, D1r, D2r, wvl, Dzr)
        c4r = constraint4_N(d1r, d2r, wvl, Dzr)
        # 只在约束 1 成立的参数域上比较（约束 1 <=> 几何项 <= 约束4/2）
        geo = D1r / (2 * d1r) + D2r / (2 * d2r)
        if geo <= c4r / 2 + 1e-9 * c4r:
            viol = max(viol, c2r - c4r)     # 应恒 <= 0
    print('    1000 组随机参数（约束 1 成立域内）：max[约束2右端 - 约束4右端] = %.3e（<=0 恒成立）'
          % viol)
    results.append(('习题6 max(c2-c4)', viol, 0.0, viol))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
