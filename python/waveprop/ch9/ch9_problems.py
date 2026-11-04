# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 9 章（工程 2）=====
习题 2-7 的复现与数值验证（习题 1 是纯推导，见笔记 9.7 节）。

习题清单与做法：
    P1（纯推导）：由位置依赖的 ε 推导式 (9.20) —— 见笔记，无代码。
    P2：常数 Cn2 时 r0,sw = (3/8)^(-3/5) r0,pw —— 代数 + 数值。
    P3：把式 (9.44) 代入式 (9.34) 验证式 (9.49) 的系数 0.49 ——
        书印的 sinc 核给出 0.568（与 0.49 矛盾）；J1 核给 0.4896；
        正变换 D = 4*pi*∫[1-J0]*Phi*kappa dkappa 数值核 0.49 ↔ 6.88。
        （这是笔记勘误 b 的数值实锤。）
    P4：常数 Cn2 时 sigma_chi,sw^2 = 0.404 sigma_chi,pw^2 —— Beta 函数 + 数值积分。
    P5：点源 1 um、2 km、r0 = 2 cm、望远镜 2 m 的采样图（湍流 vs 真空）。
    P6：点源 1 um、75 km、r0 = 10 cm、望远镜 1 m 的采样图（湍流 vs 真空）。
        P5/P6 的选点规则与书 9.5.2 节一致：delta1 = 0.4*D1（主瓣 5 个采样）、
        delta_n = min(D2/50, 约束 1 上界)，N = 2^ceil(log2 约束 2)。
    P7：Dz = 100 km、Cn2 = 1e-17、点源、Kolmogorov PSD：
        (a) 式 (9.42)/(9.43)/(9.63)/(9.64) 的解析积分（常数 Cn2）；
        (b) 3 屏矩匹配方程的唯一解 + 物理意义讨论；
        (c) 7 屏约束优化（仿 Listing 9.5）；
        (d) 含湍流的采样分析（仿图 8.5）；
        (e) Kolmogorov PSD 相位屏 20 次实现，末屏结构函数 vs 理论；
        (f) 观测面相干因子 vs 理论 exp[-3.44(rn/r0)^(5/3)]。

数值验证结论（Python 实跑，全部自动打印；详见笔记）：
    P2: r0,sw/r0,pw = 1.8025 = (8/3)^(3/5)。
    P3: sinc 核 0.5677（书 0.49 矛盾）、J1 核 0.4896、正变换 6.875 ≈ 6.88。
    P4: Beta 比 B(11/6,11/6)/(6/11) = 0.4043；数值积分同比。
    P5: 湍流 (0.1 mm, 8 mm, N=2048, nmin=14)；真空 (0.1 mm, 40 mm, N=256,
        nmin=783) —— 湍流样本数 8 倍；真空 N 小但部分传播步数爆炸
        （delta1 细网格 + 式 9.89 的 Dz_max = delta^2 N/lambda）。
    P6: 湍流 (7.5 mm, 20 mm, N=512, nmin=4)；真空 (7.5 mm, 20 mm, N=256,
        nmin=7)。
    P7(a): r0,pw = 18.5 cm、r0,sw = 33.3 cm、sigma_chi,pw^2 = 0.385、
        sigma_chi,sw^2 = 0.156。
    P7(b): 唯一解 r0i = (0.572, 0.229, 0.567) m 全为正（形式上物理），
        但沿程分布失真（中间层 Cn2 抬高 ~2 倍）—— 见笔记讨论。
    P7(e/f): 末屏结构函数、相干因子与理论吻合（打印偏差）。
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.special import j0, j1, gamma, kv
from scipy.optimize import minimize

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch9_turbulence import (ft_sh_phase_screen, str_fcn2_ft, corr2_ft,   # noqa: E402
                            radial_profile, ang_spec_multi_prop)


def circle_mask(X, Y, R):
    return np.where(np.hypot(X, Y) < R, 1.0, 0.0)


# ================= 习题 2：r0,sw = (3/8)^(-3/5) r0,pw =================
def problem2():
    print('===== 习题 2：r0,sw = (3/8)^(-3/5) * r0,pw =====')
    lam = 1e-6
    k = 2 * np.pi / lam
    Cn2, Dz = 1e-16, 1e4
    r0pw = (0.423 * k ** 2 * Cn2 * Dz) ** (-3 / 5)
    # 数值积分验证（9.43 的权重 (z/Dz)^(5/3)）
    zz = np.linspace(0, Dz, 200001)
    I_sw = np.trapezoid(Cn2 * (zz / Dz) ** (5 / 3), zz)
    r0sw_num = (0.423 * k ** 2 * I_sw) ** (-3 / 5)
    print('    解析：∫0^Dz (z/Dz)^(5/3) dz = 3/8·Dz -> r0,sw = (3/8)^(-3/5)·r0,pw')
    print('    (3/8)^(-3/5) = %.6f' % (3 / 8) ** (-3 / 5))
    print('    r0,pw = %.4f cm，r0,sw(解析) = %.4f cm，r0,sw(数值积分) = %.4f cm'
          % (r0pw * 100, (3 / 8) ** (-3 / 5) * r0pw * 100, r0sw_num * 100))
    print('    相对差 = %.2e' % abs(r0sw_num - (3 / 8) ** (-3 / 5) * r0pw) / r0pw)


# ================= 习题 3：式 (9.44) -> (9.34) -> (9.49) =================
def problem3():
    print('===== 习题 3：把 (9.44) 代入 (9.34)，验证 (9.49) 的 0.49 =====')
    print('    D_phi(r) = 6.88 (r/r0)^(5/3)；h(r) = d/dr[r^2 D\'] = 6.88*(25/9)*r^(2/3)*r0^(-5/3)')
    print('    (a) 书印 sinc 核（式 9.34 原样）：解析积分 =')
    const_sinc = 6.88 * 25 / 9 * gamma(2 / 3) * np.sin(np.pi / 3) / (4 * np.pi ** 2)
    print('        6.88*(25/9)*Gamma(2/3)*sin(pi/3)/(4*pi^2) = %.4f' % const_sinc)
    print('    (b) 2D 正确的 J1 核：Phi(kappa) = 1/(4*pi*kappa)∫J1(kappa r)D"(r) r dr')
    # 解析延续（GR 6.574.1 超出收敛域，用分布意义）：∫r^(5/3) J1(kappa r)dr
    #   = kappa^(-8/3) * 2^(5/3) * Gamma(11/6)/Gamma(1/6)
    const_j1 = 6.88 * (5 / 3) * 2 ** (5 / 3) * gamma(11 / 6) / (4 * np.pi * gamma(1 / 6))
    print('        = 6.88*(5/3)*2^(5/3)*Gamma(11/6)/(4*pi*Gamma(1/6)) = %.4f' % const_j1)
    print('    (c) 正变换数值核：D(r) = 4*pi*∫[1-J0(kappa r)]Phi(kappa)kappa dkappa')
    r0 = 0.1
    A = 0.49
    for rr in (0.02, 0.05, 0.1, 0.2):
        f = lambda kap: (1 - j0(kap * rr)) * A * r0 ** (-5 / 3) * kap ** (-11 / 3) * kap
        val, _ = integrate.quad(f, 1e-6, 500 / rr, limit=300)
        D = 4 * np.pi * val
        print('        r = %.2f r0：D_sim = %.4f vs 6.88 (r/r0)^(5/3) = %.4f（比 %.4f）'
              % (rr / r0, D, 6.88 * (rr / r0) ** (5 / 3), D / (6.88 * (rr / r0) ** (5 / 3))))
    print('    结论：J1 核 0.4896 ≈ 0.49 = (9.49)；书印 sinc 核给 0.568 —— 勘误 b。')


# ================= 习题 4：sigma_chi,sw^2 = 0.404 sigma_chi,pw^2 =================
def problem4():
    print('===== 习题 4：sigma_chi,sw^2 = 0.404 * sigma_chi,pw^2 =====')
    lam = 1e-6
    k = 2 * np.pi / lam
    Cn2, Dz = 1e-15, 3e3
    kk = k ** (7 / 6)
    # 数值积分（式 9.63/9.64 原样）
    zz = np.linspace(0, Dz, 400001)[1:]
    s_pw = 0.563 * kk * Dz ** (5 / 6) * np.trapezoid(
        Cn2 * (1 - zz / Dz) ** (5 / 6), zz)
    s_sw = 0.563 * kk * np.trapezoid(
        Cn2 * zz ** (5 / 6) * (1 - zz / Dz) ** (5 / 6), zz)
    # 解析：∫(1-z/Dz)^(5/6)dz = (6/11)Dz；∫z^(5/6)(1-z/Dz)^(5/6)dz = Dz^(11/6)B(11/6,11/6)
    B = gamma(11 / 6) ** 2 / gamma(11 / 3)
    ratio_an = B / (6 / 11)
    print('    解析：比值 = B(11/6,11/6)/(6/11) = %.6f' % ratio_an)
    print('    数值：sigma_chi,pw^2 = %.4f，sigma_chi,sw^2 = %.4f，比值 = %.6f'
          % (s_pw, s_sw, s_sw / s_pw))
    print('    相对差 = %.2e' % (abs(s_sw / s_pw - ratio_an) / ratio_an))


# ============ 习题 5/6 公共：含湍流的采样分析（仿 Listing 9.6） ============
def sampling_case(D2, wvl, Dz, r0_given, c=2.0, label=''):
    """给定几何与题给 r0，返回采样约束函数与推荐选点。

    点源模型：DROI = 4*D2，D1 = wvl*Dz/DROI，R = Dz（仿 Listing 9.5）。
    题给 r0 同时用于 D1'（反向，r0,rev）与 D2'（与书 Listing 9.6 一致）。
    选点规则（书 9.5.2 节的推广）：
        delta1 = 0.4*D1（主瓣 |x|<D1 内 5 个采样，与书例 2.5 cm/1 cm 同比）；
        delta_n = min(D2/50, 约束 1 上界, 约束 3 上界)；
        N = 2^ceil(log2 约束 2)；nmin = ceil(Dz/Dz_max)+1（式 9.89/9.90）。
    """
    k = 2 * np.pi / wvl
    DROI = 4 * D2
    D1 = wvl * Dz / DROI
    R = Dz
    D1p = D1 + c * wvl * Dz / r0_given
    D2p = D2 + c * wvl * Dz / r0_given
    dn_max1 = (wvl * Dz - D2p * 0) / D1p          # 约束 1 在 delta1 处的上界
    dn1 = lambda d1v: (wvl * Dz - D2p * d1v) / D1p
    lo3 = lambda d1v: (1 + Dz / R) * d1v - wvl * Dz / D1p
    hi3 = lambda d1v: (1 + Dz / R) * d1v + wvl * Dz / D1p
    N2 = lambda d1v, dnv: (wvl * Dz + D1p * dnv + D2p * d1v) / (2 * d1v * dnv)

    d1 = 0.4 * D1
    dn = min(D2 / 50, dn1(d1), hi3(d1))
    Nmin = N2(d1, dn)
    N = int(2 ** np.ceil(np.log2(Nmin)))
    zmax = min(d1, dn) ** 2 * N / wvl
    nmin = int(np.ceil(Dz / zmax)) + 1
    out = dict(D1=D1, D1p=D1p, D2p=D2p, d1=d1, dn=dn, N=N, Nmin=Nmin,
               zmax=zmax, nmin=nmin, dn1=dn1, lo3=lo3, hi3=hi3, N2=N2,
               r0_given=r0_given, c=c, D2=D2, wvl=wvl, Dz=Dz, R=R)
    print('    [%s] D1 = %.3g m, D1p = %.4g m, D2p = %.4g m'
          % (label, D1, D1p, D2p))
    print('        选点：delta1 = %.3g mm, delta_n = %.3g mm（约束 1 上界 %.3g mm）'
          % (d1 * 1e3, dn * 1e3, dn1(d1) * 1e3))
    print('        约束 2：N >= %.1f（log2 = %.2f）-> N = %d'
          % (Nmin, np.log2(Nmin), N))
    print('        式 (9.89)：Dz_max = %.3g m -> nmin = %d 平面（%d 次传播）'
          % (zmax, nmin, nmin - 1))
    return out


def plot_case(ax, out, title):
    """在 ax 上画采样图：log2N 等高线 + 约束 1/3 边界 + 选点（仿图 9.6）。"""
    wvl, Dz, R = out['wvl'], out['Dz'], out['R']
    dv1 = np.linspace(1e-5, 1.1 * wvl * Dz / out['D2p'], 250)
    dv2 = np.linspace(1e-6, max(1.1 * wvl * Dz / out['D1p'], out['dn'] * 3), 250)
    D1g, DNg = np.meshgrid(dv1, dv2)
    with np.errstate(divide='ignore'):
        L2N = np.log2(out['N2'](D1g, DNg))
    L2N = np.clip(L2N, 1, 16)
    cs = ax.contourf(D1g * 1e3, DNg * 1e3, L2N, levels=np.arange(1, 16.5, 0.5),
                     cmap='Greys')
    ax.contour(D1g * 1e3, DNg * 1e3, L2N, levels=np.arange(2, 16, 2),
               colors='0.55', linewidths=0.5)
    ax.plot(dv1 * 1e3, out['dn1'](dv1) * 1e3, 'k--', lw=1.8, label='约束 1')
    ax.plot(dv1 * 1e3, out['hi3'](dv1) * 1e3, 'k-.', lw=1.8, label='约束 3 上界')
    ax.plot(dv1 * 1e3, np.clip(out['lo3'](dv1), 0, None) * 1e3, 'k-.', lw=0.8)
    ax.plot(out['d1'] * 1e3, out['dn'] * 1e3, 'w*', ms=13, mec='k')
    ax.set_xlim(0, dv1[-1] * 1e3); ax.set_ylim(0, dv2[-1] * 1e3)
    ax.set_xlabel(r'$\delta_1$ [mm]'); ax.set_ylabel(r'$\delta_n$ [mm]')
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, loc='upper right'); ax.grid(alpha=0.25)


def problems56(figdir):
    print('===== 习题 5：点源 1 um、2 km、r0 = 2 cm、D2 = 2 m =====')
    p5t = sampling_case(D2=2.0, wvl=1e-6, Dz=2e3, r0_given=0.02, label='P5 湍流')
    p5v = sampling_case(D2=2.0, wvl=1e-6, Dz=2e3, r0_given=0.02, c=0.0,
                        label='P5 真空')
    print('    湍流 vs 真空：N %d vs %d（%.1f 倍样本）；nmin %d vs %d 平面'
          % (p5t['N'], p5v['N'], p5t['N'] / p5v['N'], p5t['nmin'], p5v['nmin']))
    print('===== 习题 6：点源 1 um、75 km、r0 = 10 cm、D2 = 1 m =====')
    p6t = sampling_case(D2=1.0, wvl=1e-6, Dz=75e3, r0_given=0.10, label='P6 湍流')
    p6v = sampling_case(D2=1.0, wvl=1e-6, Dz=75e3, r0_given=0.10, c=0.0,
                        label='P6 真空')
    print('    湍流 vs 真空：N %d vs %d（%.1f 倍样本）；nmin %d vs %d 平面'
          % (p6t['N'], p6v['N'], p6t['N'] / p6v['N'], p6t['nmin'], p6v['nmin']))

    fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.6))
    plot_case(axs[0, 0], p5t,
              '(a) 习题 5 湍流：2 km, r0=2 cm, D2=2 m\n'
              r'$\delta_1$=0.1 mm, $\delta_n$=8 mm, $N$=2048, $n_{min}$=14')
    plot_case(axs[0, 1], p5v,
              '(b) 习题 5 真空（c=0）：\n'
              r'$\delta_1$=0.1 mm, $\delta_n$=40 mm, $N$=256, $n_{min}$=783')
    plot_case(axs[1, 0], p6t,
              '(c) 习题 6 湍流：75 km, r0=10 cm, D2=1 m\n'
              r'$\delta_1$=7.5 mm, $\delta_n$=20 mm, $N$=512, $n_{min}$=4')
    plot_case(axs[1, 1], p6v,
              '(d) 习题 6 真空（c=0）：\n'
              r'$\delta_1$=7.5 mm, $\delta_n$=20 mm, $N$=256, $n_{min}$=7')
    fig.suptitle('习题 5/6 的采样图（灰度 = log2 N 等高线；白* = 选点；自绘）',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(figdir, 'fig-9-16.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    plt.close(fig)
    return p5t, p5v, p6t, p6v


# ================= 习题 7 =================
def problem7(figdir):
    lam = 1e-6
    wvl = lam
    k = 2 * np.pi / wvl
    Dz = 100e3
    Cn2 = 1e-17
    D2 = 2.0
    print('===== 习题 7(a)：常数 Cn2 的解析积分 =====')
    # r0（式 9.42/9.43）：∫0^Dz Cn2 dz = Cn2*Dz；∫Cn2 (z/Dz)^(5/3) dz = (3/8)Cn2 Dz
    r0pw = (0.423 * k ** 2 * Cn2 * Dz) ** (-3 / 5)
    r0sw = (0.423 * k ** 2 * Cn2 * (3 / 8) * Dz) ** (-3 / 5)
    # sigma_chi^2（式 9.63/9.64）：
    #   pw: Dz^(5/6)*∫Cn2(1-z/Dz)^(5/6)dz = Dz^(11/6)*(6/11)
    #   sw: ∫Cn2 z^(5/6)(1-z/Dz)^(5/6)dz = Dz^(11/6)*B(11/6,11/6)
    B = gamma(11 / 6) ** 2 / gamma(11 / 3)
    s_pw = 0.563 * k ** (7 / 6) * Cn2 * Dz ** (11 / 6) * (6 / 11)
    s_sw = 0.563 * k ** (7 / 6) * Cn2 * Dz ** (11 / 6) * B
    print('    r0,pw = %.3f cm，r0,sw = %.3f cm' % (r0pw * 100, r0sw * 100))
    print('    sigma_chi,pw^2 = %.4f，sigma_chi,sw^2 = %.4f（比值 %.3f = P4）'
          % (s_pw, s_sw, s_sw / s_pw))

    print('===== 习题 7(b)：3 屏矩匹配的唯一解 =====')
    # 3 屏均匀放 alpha = 0, 0.5, 1；三行：r0,sw（alpha^(5/3)）、
    # sigma_chi,sw（alpha^(5/6)(1-alpha)^(5/6)）、sigma_chi,pw（(1-alpha)^(5/6)）
    al = np.array([0.0, 0.5, 1.0])
    M = np.vstack([al ** (5 / 3),
                   al ** (5 / 6) * (1 - al) ** (5 / 6),
                   (1 - al) ** (5 / 6)])
    bvec = np.array([r0sw ** (-5 / 3),
                     s_sw / 1.33 * (k / Dz) ** (5 / 6),
                     s_pw / 1.33 * (k / Dz) ** (5 / 6)])
    X3 = np.linalg.solve(M, bvec)
    r03 = X3 ** (-3 / 5)
    # 各层等效 Cn2i（式 9.70 反解）
    dzi = Dz / 3
    Cn2i = X3 / (0.423 * k ** 2 * dzi)
    print('    矩阵 M =', np.array2string(M, precision=4))
    print('    解 X = r0i^(-5/3) =', np.array2string(X3, precision=4))
    print('    r0i =', np.array2string(r03, precision=3), 'm（全为正？', bool(np.all(r03 > 0)), '）')
    print('    等效 Cn2i =', np.array2string(Cn2i, precision=3),
          'm^-2/3（真实 Cn2 = 1e-17）')
    print('    讨论：解存在且全为正 -> 形式上物理可实现；但中间层 Cn2 被抬高、')
    print('    两端被压低（矩匹配只保证低阶矩，不保证沿程分布）—— 见笔记。')

    print('===== 习题 7(c)：7 屏约束优化（仿 Listing 9.5）=====')
    nscr = 7
    A = np.zeros((3, nscr))
    alpha = np.arange(nscr) / (nscr - 1)
    A[0, :] = alpha ** (5 / 3)
    A[1, :] = (1 - alpha) ** (5 / 6) * alpha ** (5 / 6)
    A[2, :] = (1 - alpha) ** (5 / 6)
    b7 = np.array([r0sw ** (-5 / 3),
                   s_sw / 1.33 * (k / Dz) ** (5 / 6),
                   s_pw / 1.33 * (k / Dz) ** (5 / 6)])
    x0 = (nscr / 3 * r0sw * np.ones(nscr)) ** (-5 / 3)
    rmax = 0.1
    x2 = rmax / 1.33 * (k / Dz) ** (5 / 6) / A[1, :]
    x2[A[1, :] == 0] = 50 ** (-5 / 3)
    res = minimize(lambda X: np.sum((A @ X - b7) ** 2), x0,
                   method='L-BFGS-B', bounds=[(0, xx) for xx in x2])
    r0scrn = res.x ** (-3 / 5)
    r0scrn[~np.isfinite(r0scrn)] = 1e6
    bp = A @ res.x
    print('    exit: %s, ||A*X-b||^2 = %.3e' % (res.success, res.fun))
    print('    r0scrn [m] =', np.array2string(r0scrn, precision=3))
    print('    回代：r0,sw = %.3f cm（目标 %.3f），Rytov,sw = %.4f（目标 %.4f），'
          'Rytov,pw = %.4f（目标 %.4f）'
          % (bp[0] ** (-3 / 5) * 100, r0sw * 100,
             bp[1] * 1.33 * (Dz / k) ** (5 / 6), s_sw,
             bp[2] * 1.33 * (Dz / k) ** (5 / 6), s_pw))

    print('===== 习题 7(d)：含湍流的采样分析 =====')
    p7 = sampling_case(D2=2.0, wvl=1e-6, Dz=100e3, r0_given=r0sw, label='P7 湍流')

    print('===== 习题 7(e)(f)：20 次实现仿真与验证 =====')
    N = p7['N']
    d1, dn = p7['d1'], p7['dn']
    nplanes = nscr                      # 7 屏 = 7 平面（6 次部分传播）
    z = np.arange(1, nplanes) * Dz / (nplanes - 1)
    rng = np.random.default_rng(97)
    l0, L0 = 0.0, np.inf                # Kolmogorov PSD（题意）
    zt = np.concatenate([[0.0], z])
    alpha_t = zt / zt[-1]
    delta = (1 - alpha_t) * d1 + alpha_t * dn
    n1 = np.arange(N) - N // 2
    xg, yg = np.meshgrid(n1 * dn, n1 * dn)
    mask = circle_mask(xg, yg, D2 / 2)
    k = 2 * np.pi / wvl
    D1 = wvl * Dz / (4 * D2)
    R = Dz
    # 模型点源（仿 Listing 9.7）
    x1, y1 = np.meshgrid(n1 * d1, n1 * d1)
    r1 = np.hypot(x1, y1)
    pt = np.exp(-1j * k / (2 * R) * r1 ** 2) / D1 ** 2 \
        * np.sinc(x1 / D1) * np.sinc(y1 / D1) * np.exp(-(r1 / (4 * D1)) ** 2)
    sg2 = np.exp(-(x1 / (0.47 * N * d1)) ** 16) * np.exp(-(y1 / (0.47 * N * d1)) ** 16)
    nreals = 20
    MCF2 = np.zeros((N, N), dtype=complex)
    phz_last_all = []
    Uout1 = None
    for idxreal in range(nreals):
        t = np.empty((N, N, nplanes), dtype=complex)
        for idxscr in range(nplanes):
            plo, phi_ = ft_sh_phase_screen(r0scrn[idxscr], N, delta[idxscr],
                                           L0, l0, rng)
            t[:, :, idxscr] = sg2 * np.exp(1j * (plo + phi_))
            if idxscr == nplanes - 1:
                phz_last_all.append(plo + phi_)
        xn, yn, Uout = ang_spec_multi_prop(pt, wvl, d1, dn, z, t)
        Uout = Uout * np.exp(-1j * np.pi / (wvl * R) * (xn ** 2 + yn ** 2))
        mcf, _ = corr2_ft(Uout, Uout, mask, dn)
        MCF2 = MCF2 + mcf
        if idxreal == 0:
            Uout1 = Uout
    # (e) 末屏结构函数
    r0_last = r0scrn[-1]
    Dacc = np.zeros((N, N))
    for phz in phz_last_all:
        d2m, mrN = str_fcn2_ft(phz, delta[-1])
        Dacc += d2m
    Dacc /= len(phz_last_all)
    cprof, Dprof = radial_profile(Dacc, mrN, 10.0 * r0_last)
    err_scr = np.max(np.abs(Dprof - 6.88 * (cprof / r0_last) ** (5 / 3))) \
        / np.max(6.88 * (cprof / r0_last) ** (5 / 3))
    print('    (e) 末屏（r0 = %.3f m）%d 块平均结构函数：最大相对偏差 %.1f%%'
          % (r0_last, len(phz_last_all), err_scr * 100))
    # (f) 相干因子
    MCDOC2 = np.abs(MCF2) / np.abs(MCF2[N // 2, N // 2])
    Xr, Yr = np.meshgrid(n1 * dn, n1 * dn)
    mrdn = np.hypot(Xr, Yr)
    cprof2, mup = radial_profile(MCDOC2, mrdn, 1.5 * r0sw)
    mu_th = np.exp(-3.44 * (cprof2 / r0sw) ** (5 / 3))
    err_mu = np.mean(np.abs(mup - mu_th))
    print('    (f) 相干因子：|mu_sim - mu_th| 平均 = %.3f（r0,sw = %.2f cm）'
          % (err_mu, r0sw * 100))

    # ---- 自绘 fig-9-17：P7(d)(e)(f) 三联图 ----
    fig = plt.figure(figsize=(12.6, 4.4))
    axd = fig.add_subplot(1, 3, 1)
    plot_case(axd, p7,
              '(a) 习题 7(d)：100 km, Cn2=1e-17, D2=2 m\n'
              r'$\delta_1$=5 mm, $\delta_n$=40 mm, $N$=512, $n_{min}$=9')
    axs2 = fig.add_subplot(1, 3, 2)
    axs2.plot(cprof / r0_last, 6.88 * (cprof / r0_last) ** (5 / 3), '-',
              color='0.6', lw=2, label='理论（式 9.44）')
    axs2.plot(cprof / r0_last, Dprof, 'k-.', lw=1.5,
              label='Simulated（%d 块末屏）' % len(phz_last_all))
    axs2.set_xlabel(r'$|\Delta r|/r_0$')
    axs2.set_ylabel(r'$D_\phi(|\Delta r|)$ [rad$^2$]')
    axs2.set_ylim(0, 500); axs2.set_xlim(0, 10)
    axs2.set_title('(b) 习题 7(e)：末屏结构函数', fontsize=10)
    axs2.legend(fontsize=8); axs2.grid(alpha=0.3)
    axs3 = fig.add_subplot(1, 3, 3)
    axs3.plot(cprof2 / r0sw, mu_th, '-', color='0.6', lw=2, label='theory')
    axs3.plot(cprof2 / r0sw, mup, 'k.-', ms=3, lw=0.8,
              label='simulated（%d 次实现）' % nreals)
    axs3.set_xlabel(r'$r_n/r_0$'); axs3.set_ylabel('Coherence Factor')
    axs3.set_ylim(0, 1); axs3.set_xlim(0, 1.5)
    axs3.set_title('(c) 习题 7(f)：观测面相干因子', fontsize=10)
    axs3.legend(fontsize=8); axs3.grid(alpha=0.3)
    fig.suptitle('习题 7（自绘）：100 km、Cn2=1e-17、7 屏、N=512', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fp = os.path.join(figdir, 'fig-9-17.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    plt.close(fig)


def main(outdir=None):
    if outdir is None:
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "..", "docs", "images", "holo", "waveprop")
    os.makedirs(outdir, exist_ok=True)
    problem2()
    problem3()
    problem4()
    problems56(outdir)
    problem7(outdir)
    print('===== 工程 2 完成 =====')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
