# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 9 章（工程 2）：习题复现 =====
习题 1-7 可数值化部分的复现（原书 9.7 节）；与 MATLAB 版
matlab/waveprop/ch9/ch9_problems.m 数值行为一致。

习题清单与处理：
    1. 位置相关介电常数的 Maxwell 方程组合 -> 式 (9.20)：解析推导
       （见笔记正文），本脚本不做数值。
    2. 常数 Cn2 时 r0,sw = (3/8)^(-3/5) r0,pw：数值验证。
    3. 把式 (9.44) 代入式 (9.34) 检验相位 PSD 系数：数值积分。
       书印核 sin(kr)/(kr)（3D 场公式）与相位这种 2D 随机场的
       逆变换核 J1 不一致——J1 核给 0.49（= 式 9.49），sin 核给
       ~0.35（勘误 b 实锤）。Kolmogorov 的 r^(5/3) 结构函数使两个
       积分都只能在阻尼/Abel 求和意义下取值：统一用
       I(eps) = ∫...e^(-eps r)dr，eps 减半外推 eps->0。
    4. 常数 Cn2 时 sigma_chi^2,sw = 0.404 sigma_chi^2,pw：数值验证。
    5. 点源 lambda=1um、2 km、r0=2 cm、望远镜 2 m：真空 vs 湍流
       采样图（Coy 展宽模型），比较网格点数与部分传播段数。
    6. 同 5：75 km、r0=10 cm、望远镜 1 m。
    7. 点源 lambda=1um、Dz=100 km、Cn2=1e-17（Kolmogorov 谱）：
       (a) 解析算 r0 与 log 振幅方差（平面波/点源，常数 Cn2）；
       (b) 3 屏 3 方程的唯一解及其物理性讨论；
       (c) 7 屏的有界最小二乘（对应 Listing 9.5 的 fmincon 流程）；
       (d) 采样图（接收口径按题 2 m；为控制算力，传播验证的
           "示例规模" 用 D2 = 0.5 m、N=512，见 ch9_turbulence.py p7）；
       (e)(f) 与 9.5 节同流程的相位屏/相干因子验证，由
           `python ch9_turbulence.py p7` 完成（20 次实现）。

用法：python ch9_problems.py [outdir]   （输出 out_prob*.png 到 outdir）
依赖：numpy、matplotlib；ch9_turbulence（同目录）。
"""

import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import j1 as besselj1

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'dejavusans'   # 数学字体用 DejaVu（含 U+2212）

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch9_turbulence import atmos_setup  # noqa: E402


def prob2():
    """习题 2：r0,sw/r0,pw = (3/8)^(-3/5)。"""
    wvl, Dz, Cn2 = 1e-6, 50e3, 1e-16
    k = 2 * np.pi / wvl
    r0pw = (0.423 * k ** 2 * Cn2 * Dz) ** (-3.0 / 5)
    r0sw = (0.423 * k ** 2 * Cn2 * 3 / 8 * Dz) ** (-3.0 / 5)
    th = (3 / 8) ** (-3.0 / 5)
    print('===== 习题 2：r0,sw/r0,pw =====')
    print('    r0,pw = %.4f m, r0,sw = %.4f m；比值 = %.6f，'
          '理论 (3/8)^(-3/5) = %.6f，差 %.1e'
          % (r0pw, r0sw, r0sw / r0pw, th, abs(r0sw / r0pw - th)))
    return r0sw / r0pw


def _damped_int(fun, eps, rmax, dr):
    r = np.arange(dr, rmax, dr)
    return np.trapezoid(fun(r) * np.exp(-eps * r), dx=dr)


def prob3():
    """习题 3：D_K 代入式 (9.34) 检验相位 PSD 系数（勘误 b）。
    r0 = 1 m、D(r) = 6.88 r^(5/3)；在 kappa = 1 rad/m 处算 Phi(1)。
    两个积分都按 I(eps)（eps 减半）外推 eps -> 0（Abel 求和）。"""
    dp = lambda r: 6.88 * (5.0 / 3) * r ** (2.0 / 3)          # D'(r)
    dr2 = lambda r: 6.88 * (5.0 / 3) * (8.0 / 3) * r ** (5.0 / 3)  # d/dr(r^2 D')
    kap = 1.0
    epss = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    # 2D 正确核（J1）：Phi = 1/(4 pi kappa) ∫ J1(kappa r) D'(r) r dr
    vJ = [_damped_int(lambda r: besselj1(kap * r) * dp(r) * r, e, 6000.0, 2e-3)
          for e in epss]
    pJ = np.polyfit(epss, vJ, 1)                    # 线性外推 eps->0
    Phi_J1 = pJ[1] / (4 * np.pi * kap)
    # 解析（Abel）：∫ u^(5/3) J1(u) du = 2^(5/3) Gamma(11/6)/Gamma(1/6)
    J1_exact = 2 ** (5.0 / 3) * math.gamma(11.0 / 6) / math.gamma(1.0 / 6)
    Phi_J1_exact = 6.88 * (5.0 / 3) * J1_exact / (4 * np.pi) * kap ** (-8.0 / 3)
    # 书印 3D 核（sin(kr)/(kr)）：Phi = 1/(4 pi^2 kappa^2) ∫ sin(kr)/(kr) d/dr(r^2 D') dr
    vS = [_damped_int(lambda r: np.sin(kap * r) / (kap * r) * dr2(r),
                      e, 6000.0, 2e-3) for e in epss]
    pS = np.polyfit(epss, vS, 1)
    Phi_sin = pS[1] / (4 * np.pi ** 2 * kap ** 2)
    S_exact = math.gamma(5.0 / 3) * math.sin(5 * np.pi / 6)   # ∫ r^(2/3) sin r dr
    Phi_sin_exact = (6.88 * (5.0 / 3) * (8.0 / 3) * S_exact
                     / (4 * np.pi ** 2)) * kap ** (-11.0 / 3)
    print('===== 习题 3：式 (9.44) -> 式 (9.34) 的相位 PSD 系数 =====')
    print('    J1 核（2D 正确逆变换）：Phi(k=1) 数值外推 = %.4f，解析（Abel）= %.4f'
          ' —— 吻合式 (9.49) 的 0.49' % (Phi_J1, Phi_J1_exact))
    print('    sin 核（书印 3D 核）：数值外推 = %.4f，解析 = %.4f —— 不是 0.49'
          '（式 9.34 的核应为 J1，勘误 b）' % (Phi_sin, Phi_sin_exact))
    return Phi_J1, Phi_sin


def prob4():
    """习题 4：sigma_chi^2,sw / sigma_chi^2,pw = 0.404。"""
    # 常数 Cn2：式 (9.63) 积分 = (6/11) L^(11/6)；式 (9.64) = B(11/6,11/6) L^(11/6)
    B = math.exp(2 * math.lgamma(11.0 / 6) - math.lgamma(11.0 / 3))
    ratio = B / (6.0 / 11)
    print('===== 习题 4：sigma_sw^2 / sigma_pw^2 =====')
    print('    B(11/6, 11/6) = %.6f；比值 = %.6f（理论 0.404），差 %.1e'
          % (B, ratio, abs(ratio - 0.404)))
    return ratio


def _samp_case(Dz, r0, D2, lam=1e-6):
    """习题 5/6：真空 vs 湍流的采样账本（Coy 展宽模型 c=2）。"""
    DROI = 4 * D2
    D1 = lam * Dz / DROI              # 点源模型中央瓣全宽 2*D1
    delta1 = 2 * D1 / 5               # 中央瓣 5 点采样（书例做法）
    # ---- 真空（第 8 章约束）----
    dn_v = min((lam * Dz - D2 * delta1) / D1, D2 / 25)    # 约束 1 + 口径 25 点
    N2_v = (D1 / (2 * delta1) + D2 / (2 * dn_v)
            + lam * Dz / (2 * delta1 * dn_v))             # 约束 2
    N_v = 2 ** np.ceil(np.log2(N2_v))
    zmax_v = min(delta1, dn_v) ** 2 * N_v / lam           # 式 (9.89)
    n_v = int(np.ceil(Dz / zmax_v)) + 1                   # 式 (9.90)
    # ---- 湍流（式 9.84-9.88，c=2）----
    c = 2
    D1p = D1 + c * lam * Dz / r0
    D2p = D2 + c * lam * Dz / r0
    c1_t = (lam * Dz - D2p * delta1) / D1p                # 约束 1'（式 9.86）
    dn_t = min(c1_t, D2 / 25)
    N2_t = (D1p / (2 * delta1) + D2p / (2 * dn_t)
            + lam * Dz / (2 * delta1 * dn_t))             # 约束 2'（式 9.87）
    N_t = 2 ** np.ceil(np.log2(N2_t))
    zmax_t = min(delta1, dn_t) ** 2 * N_t / lam
    n_t = int(np.ceil(Dz / zmax_t)) + 1
    return dict(D1=D1, delta1=delta1, dn_v=dn_v, N2_v=N2_v, N_v=N_v,
                zmax_v=zmax_v, n_v=n_v, D1p=D1p, D2p=D2p, c1_t=c1_t, dn_t=dn_t,
                N2_t=N2_t, N_t=N_t, zmax_t=zmax_t, n_t=n_t)


def prob56(outdir):
    """习题 5/6：采样图与账本。"""
    os.makedirs(outdir, exist_ok=True)
    cases = [('习题 5：λ=1 μm, Δz=2 km, r0=2 cm, D2=2 m', 2e3, 0.02, 2.0),
             ('习题 6：λ=1 μm, Δz=75 km, r0=10 cm, D2=1 m', 75e3, 0.10, 1.0)]
    print('===== 习题 5/6：真空 vs 湍流的采样分析（c=2）=====')
    res = []
    for lab, Dz, r0, D2 in cases:
        s = _samp_case(Dz, r0, D2)
        res.append((lab, s))
        print('    %s' % lab)
        print('      点源模型：D1 = %.4f m，delta1 = %.4f m（瓣内 5 点）'
              % (s['D1'], s['delta1']))
        print('      真空：deltan = %.4f m；N >= %.1f -> N = 2^%d = %d；'
              'Dz_max = %.2f km -> 最少 %d 平面'
              % (s['dn_v'], s['N2_v'], np.log2(s['N_v']), s['N_v'],
                 s['zmax_v'] / 1e3, s['n_v']))
        print("      湍流：D1' = %.3f m, D2' = %.3f m；约束 1' 上界 %.4f m；"
              'N >= %.1f -> N = 2^%d = %d；Dz_max = %.2f km -> 最少 %d 平面'
              % (s['D1p'], s['D2p'], s['c1_t'], s['N2_t'], np.log2(s['N_t']),
                 s['N_t'], s['zmax_t'] / 1e3, s['n_t']))
    # 采样图：log2 N 等高线 + 约束 1 上界 + 选中点（真空/湍流两版）
    lam = 1e-6
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.0))
    for ax, (lab, Dz, r0, D2), s in zip(axs, cases, [r[1] for r in res]):
        D1p, D2p = s['D1p'], s['D2p']
        d1 = np.linspace(0.05 * s['delta1'], 1.1 * lam * Dz / D2p, 300)
        dn = np.linspace(0.05 * s['dn_t'], 1.1 * lam * Dz / D1p, 300)
        D1g, D2g = np.meshgrid(d1, dn)
        N2g = (lam * Dz + D1p * D2g + D2p * D1g) / (2 * D1g * D2g)
        cs = ax.contourf(D1g * 1e3, D2g * 1e3, np.log2(np.maximum(N2g, 1)),
                         levels=np.arange(4, 16, 1), cmap='Greys')
        fig.colorbar(cs, ax=ax, label='log2 N')
        ax.plot(d1 * 1e3, (lam * Dz - D2p * d1) / D1p * 1e3, 'k--', lw=2,
                label="约束 1'（湍流）")
        ax.plot(s['delta1'] * 1e3, s['dn_t'] * 1e3, 'r*', ms=14, zorder=5)
        ax.set_xlim(0, d1.max() * 1e3)
        ax.set_ylim(0, dn.max() * 1e3)
        ax.set_xlabel(r'$\delta_1$ [mm]')
        ax.set_ylabel(r'$\delta_n$ [mm]')
        ax.set_title(lab + '\n红* = 湍流账本选中点，N=2^%d' % np.log2(s['N_t']),
                     fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
    fig.suptitle('图 9-16（自绘仿真）习题 5/6 的采样图：湍流展宽抬高约束 1 上界，'
                 "但 D1' 变大要求更多网格点", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fp = os.path.join(outdir, 'out_prob56_sampling.png')
    fig.savefig(fp, dpi=150)
    plt.close(fig)
    print('    已保存 %s' % fp)


def prob7abcd(outdir):
    """习题 7(a)-(d)。"""
    os.makedirs(outdir, exist_ok=True)
    lam, Dz, Cn2 = 1e-6, 100e3, 1e-17
    print('===== 习题 7(a)：解析 r0 与 log 振幅方差 =====')
    st = atmos_setup(D2=0.5, wvl=lam, Dz=Dz, Cn2=Cn2, nscr=7)
    k = 2 * np.pi / lam
    r0pw = (0.423 * k ** 2 * Cn2 * Dz) ** (-3.0 / 5)
    r0sw = (0.423 * k ** 2 * Cn2 * 3 / 8 * Dz) ** (-3.0 / 5)
    B = math.exp(2 * math.lgamma(11.0 / 6) - math.lgamma(11.0 / 3))
    s_pw = 0.563 * k ** (7.0 / 6) * Cn2 * Dz ** (11.0 / 6) * (6.0 / 11)
    s_sw = 0.563 * k ** (7.0 / 6) * Cn2 * Dz ** (11.0 / 6) * B
    print('    r0,pw = %.4f m, r0,sw = %.4f m（比值 %.4f = (3/8)^(-3/5)）'
          % (r0pw, r0sw, r0sw / r0pw))
    print('    sigma_chi^2,pw = %.4f（Rytov 数 4*sigma^2 = %.4f）、'
          'sigma_chi^2,sw = %.4f' % (s_pw, 4 * s_pw, s_sw))

    print('===== 习题 7(b)：3 屏 3 方程的唯一解 =====')
    a3 = np.array([0.0, 0.5, 1.0])
    A = np.vstack([a3 ** (5.0 / 3),
                   (1 - a3) ** (5.0 / 6) * a3 ** (5.0 / 6),
                   (1 - a3) ** (5.0 / 6)])
    b = np.array([r0sw ** (-5.0 / 3),
                  s_sw / 1.33 * (k / Dz) ** (5.0 / 6),
                  s_pw / 1.33 * (k / Dz) ** (5.0 / 6)])
    X = np.linalg.solve(A, b)
    r0_3 = X ** (-3.0 / 5)
    print('    方程（行 = r0,sw / Rytov,sw / Rytov,pw；列 = 屏 1,2,3）：')
    for row, bb in zip(A, b):
        print('      [%s] * X = %.4f' % (' '.join('%7.4f' % v for v in row), bb))
    print('    唯一解 X = r0i^(-5/3) = [%s]' % ' '.join('%.3f' % v for v in X))
    print('    r0i = [%s] m' % ' '.join('%.3f' % v for v in r0_3))
    print('    物理性讨论：r0i 本身非负，但被匹配的 sigma_chi^2,pw = %.4f'
          '（Rytov 数 4*sigma^2 = %.2f >> 1）已深处强起伏区，弱起伏理论'
          '自相矛盾；且 3 屏无自由度，两端屏各被单个方程钉死，对建模误差'
          '零鲁棒。工程上应改用 7 屏 + 约束优化（见 (c)）。' % (s_pw, 4 * s_pw))

    print('===== 习题 7(c)：7 屏有界最小二乘（对应 Listing 9.5）=====')
    r0_7 = st['r0scrn']
    print('    7 屏 r0i = [%s] m' % ' '.join('%.3f' % v for v in r0_7))
    print('    回代：r0,sw = %.4f m（目标 %.4f），sigma^2 = %.4f（目标 %.4f）'
          % (st['chk'][0], r0sw, st['chk'][1], s_sw))

    print('===== 习题 7(d)：采样分析（口径 2 m 按题意）=====')
    D2p = 2.0
    DROI = 4 * D2p
    D1 = lam * Dz / DROI
    delta1 = 2 * D1 / 5
    D1p = D1 + 2 * lam * Dz / r0sw
    D2p2 = D2p + 2 * lam * Dz / r0sw
    dn_max = (lam * Dz - D2p2 * delta1) / D1p
    dn = min(dn_max, D2p / 25)
    N2 = (D1p / (2 * delta1) + D2p2 / (2 * dn) + lam * Dz / (2 * delta1 * dn))
    N = 2 ** np.ceil(np.log2(N2))
    zmax = min(delta1, dn) ** 2 * N / lam
    nmin = int(np.ceil(Dz / zmax)) + 1
    print('    D1 = %.4f m, delta1 = %.4f m；D1\' = %.3f m, D2\' = %.3f m'
          % (D1, delta1, D1p, D2p2))
    print("    约束 1'：delta_n <= %.4f m；约束 2'：N >= %.0f -> N = 2^%d = %d；"
          'Dz_max = %.1f km -> 最少 %d 平面'
          % (dn_max, N2, np.log2(N), N, zmax / 1e3, nmin))
    print('    （口径 2 m 需 N=2^%d；笔记的传播验证按示例规模取 D2 = 0.5 m'
          ' -> N=512，见 ch9_turbulence.py p7）' % np.log2(N))

    # 图：3 屏解 vs 7 屏解
    fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.4))
    axs[0].bar(np.arange(3) + 1, r0_3, color='0.6')
    axs[0].set_xticks(np.arange(3) + 1)
    axs[0].set_xlabel('相位屏编号（α = 0, 0.5, 1）')
    axs[0].set_ylabel(r'$r_{0i}$ [m]')
    axs[0].set_title('(b) 3 屏 3 方程唯一解（r0i = %.2f/%.2f/%.2f m）\n'
                     '匹配的 $\\sigma^2_{\\chi,pw}$ 对应 Rytov 数 %.2f，弱起伏理论失效'
                     % (r0_3[0], r0_3[1], r0_3[2], 4 * s_pw), fontsize=10)
    axs[1].bar(np.arange(7) + 1, np.minimum(r0_7, 2.0), color='0.6')
    axs[1].set_xticks(np.arange(7) + 1)
    axs[1].set_xlabel('相位屏编号（α = 0, 1/6, ..., 1）')
    axs[1].set_ylabel(r'$r_{0i}$ [m]')
    axs[1].set_title('(c) 7 屏有界最小二乘解（Listing 9.5 流程）\n'
                     r'回代 $r_{0,sw}/\sigma_\chi^2$ 与目标一致（两端屏权重为 0）',
                     fontsize=10)
    fig.suptitle('图 9-17（自绘仿真）习题 7(b)(c)：分层相位屏强度的两种解法',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fp = os.path.join(outdir, 'out_prob7_screens.png')
    fig.savefig(fp, dpi=150)
    plt.close(fig)
    print('    已保存 %s' % fp)


def main(argv):
    outdir = argv[1] if len(argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    prob2()
    prob3()
    prob4()
    prob56(outdir)
    prob7abcd(outdir)


if __name__ == '__main__':
    main(sys.argv)
