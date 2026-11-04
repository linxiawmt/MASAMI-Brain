# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 8 章（工程 1）=====
部分传播与放宽采样约束（Relaxed Sampling Constraints with Partial Propagations）
（原书 8.1-8.4 节，Listings 8.1-8.2）

原书要点：
    第 7 章的四条采样约束里，约束 1/2/3 由传播几何决定（Δz、D1、D2 固定时
    没得选），通常唯一的自由参数是 N；而约束 4（角谱法 N >= λΔz/(δ1δ2)）
    只依赖传播方法本身——Δz 大时它把 N 逼到离谱（N > 4096）。
    解法：把一次全长传播拆成 n-1 次"部分传播"（n 个平面），每次传播距离
    Δzi = Δz/(n-1) 变短，约束 4 按比例放松；代价是 FFT 次数变多，但换来
    N 下降，总时间通常反而更短。
    关键配套：吸收边界（8.1 节）。中间平面上把网格边缘的场乘一个衰减因子
    （超高斯 / Tukey 窗），抑制沿整个传播路径的卷绕（wrap-around）。
    网格间距沿 z 线性插值（式 8.8：δi = (1-αi)δ1 + αiδn），这使相邻两次
    部分传播在中间平面留下的两个二次相位恰好互为共轭、完全相消
    （式 8.10-8.15）——这是整套方法能自洽的代数核心。
    约束 1/2 原样保留；约束 3 代换 δ2、Δz1 后与单次传播的式 (7.53) 完全
    相同（式 8.20-8.22）；约束 4 改为对每段的要求（式 8.23），用保守界限
    Δzi <= min(δ1,δn)^2 N/λ（式 8.24）定出最少段数 n-1 = ceil(Δz/Δzmax)。

原书勘误/不一致（本脚本数值实锤，详见笔记）：
    (a) Listing 8.2 与正文/图 8.5 参数不一致：正文说"30 点跨 D2 给
        δn <= 133 µm"（即 D2 = 4 mm），图 8.5 白叉也落在 (66.7, 133) µm；
        但 Listing 8.2 印的是 D2 = 6e-3 m（δn = 200 µm）。
        图 8.6 的采样点序列（2.56/2.27/1.55/0.83/...）与 δn = 200 µm 吻合
        ——即Listing 是实际运行版本，图 8.5 白叉是旧稿残留。两组参数都可行。
    (b) 图 8.5 的约束 1 虚线截距实测 ≈ 500 µm = λΔz/(2D1)，而正文式 (7.14)
        在 Δz = 2 m、D2 = 4 mm 下应为 1000 µm（同图的约束 2 等高线
        却与 Δz = 2 m、D2 = 4 mm 精确吻合）——虚线疑为旧稿（Δz = 1 m）
        或绘图系数错误。
    (c) Listing 8.1 的吸收边界 sg 在循环里对每次部分传播（包括最后一段）
        之场都乘一遍——观测面也被衰减；式 (8.9)/(8.16)（两次传播）只在
        中间平面加 A[r2]。感兴趣区域内 sg ≈ 1 无影响，但严格做法是只在
        中间平面衰减（本脚本 absorb='between'）。
    (d) 式 (8.2) Tukey 窗的印刷 L/N 混排（条件 "r >= αL/2" 应为
        "r <= αL/2"，第二支条件/变量里 N 应为 L），且径向窗在 r > L/2
        （网格角落）需截为零。Listing 8.1 实际用的是超高斯，不受影响。
    (e) "至少 n = ceil(Δz/Δzmax)+1 次部分传播"：严格最小是
        ceil(Δz/Δzmax) 段（每段 Δz/n <= Δzmax），+1 是安全裕量——
        书例 4 段（0.5 m）即可满足，书中用了 5 段（0.4 m）。

Python 实跑结果（MATLAB 版按相同参数复现，见 matlab/waveprop/ch8/）：
    验证 0：n=2 部分传播（不加吸收边界）与第 6 章 ang_spec_prop 逐位一致
            （最大差 0，算符链自洽）。
    验证 1：书例约束总账：δ1=66.7 µm、δn=200/133 µm、N=128；
            约束 2：Listing 变体 N>=105（满足）、正文变体 N>=142.5（128
            轻微违反）；式 (8.25) Δzmax = 0.5689 m（书印 0.567 m），
            n = ceil(2/0.5689)+1 = 5 段（严格最小 4 段）。
    验证 2：Listing 8.2 复现（D2=6 mm, δn=200 µm, N=128, 5 段）：
            感兴趣区域 |dI|/峰值 = 3.9e-3、相位残差 4.8e-3 rad；
            采样点序列 I(0/0.2/0.4/0.6/0.8 mm) = 2.55/2.26/1.55/0.83/0.40
            与原书图 8.6 读数一致（也实锤 δn=200 µm 是实际运行参数）。
            正文变体（D2=4 mm, δn=133 µm）：3.9e-3 / 6.0e-3——两组都行。
    验证 3：固定网格加大 Δz（2/4/6/8 m，N=128 不变）：
            单次 |dI|/pk = 4.2e-3 → 1.1e-1（相位残差到 0.17 rad），
            多段（段距<=0.4 m）始终 <= 8.2e-3——单次退化 14 倍、多段稳定。
    验证 4：约束 3 与部分传播无关（式 8.20->8.22）：200 组随机参数，界折算
            到 δn 空间与式 (8.22) 最大差 8e-15（以区间半宽为单位），判定
            200/200 一致。补充推导：后续段（i>1、R 有限）约束中心偏移
            Δz(δi-δ1)/R，R=∞ 时严格不变（书例 R=∞ 不受影响）。
    验证 5：吸收边界勘误 (c) 实锤：book 模式（观测面也被衰减）与 between
            模式在感兴趣区域一致（sg≈1）；网格边缘 r=60 样本处振幅比
            |U_book|/|U_between| = 0.384 = sg 值。
    验证 6：Listing 8.1 的 sg = exp(-(r/0.47N)^16) 与式 (8.1) 超高斯一致
            （逐位，r^16 与 (r/σ)^16 排列差 1 ULP）；边界值：中心 1.0、
            r=0.47N 处 0.384、网格角 6e-300。

MATLAB/Python 差异（对应 .m 内注释）：
    - MATLAB 1-based / numpy 0-based：居中排列零坐标样本位置不同；
    - MATLAB 的 nsq.^8 是 double 运算；numpy 整数网格 nsq**8 会 int64
      溢出（本脚本用 float 网格规避）——复现本书的第一个坑；
    - meshgrid 行为一致；Python 复用第 6 章的 ift2 与第 7 章的约束函数；
    - 本机无 MATLAB 时以 Python 实跑结果为准。

依赖：numpy + matplotlib + scipy；fourier/common/base.py（ft2_centered）、
      waveprop/ch6（ift2/ang_spec_prop）、waveprop/ch7（约束函数）。
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ('..', '../..', '../../fourier/common', '../ch6', '../ch7'):
    _ap = os.path.abspath(os.path.join(_HERE, _p))
    if _ap not in sys.path:
        sys.path.insert(0, _ap)

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1

from base import ft2_centered as ft2, rect      # 定标正确的居中二维 DFT / rect
from ch6_fresnel_prop import ift2, ang_spec_prop, fresnel_prop_square_ap
from ch7_sampling import (constraint1_d2, constraint2_N,
                          constraint3_d2, constraint4_N)

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


# ===== 书 Listing 8.1：任意次部分传播的角谱法（真空） =====
def sg_absorbing_boundary(N, frac=0.47, n_order=8):
    """书 Listing 8.1 第 9-12 行的超高斯吸收边界（式 8.1 的实现版）。

    sg = exp(-(r/w)^16)，w = 0.47*N（以"样本数"为单位，r 为索引坐标的
    径向距离）。书代码写法 exp(-nsq.^8/w^16)（nsq = r^2，(r^2)^8 = r^16）。
    每个中间平面上网格边缘（r >= 0.47N）的场被强烈衰减，中心（ROI）处
    sg ≈ 1（例如 N=128、r=15 个样本时 sg ≈ e^-22）。
    n_order=8 对应书代码里的指数 8（即 (r/w)^16，式 8.1 的 n=16）。
    """
    n = (np.arange(N) - N // 2).astype(float)     # MATLAB: -N/2 : N/2-1（0-based 对应）
    nx, ny = np.meshgrid(n, n)                    # 与 MATLAB meshgrid 一致：nx 沿列变
    nsq = nx ** 2 + ny ** 2                       # r^2（样本单位）—— 必须是 float！
    # MATLAB 的 nsq.^8 是 double 运算（8192^8 ≈ 2e31 无压力）；numpy 的整数
    # 网格 nsq**8 会 int64 溢出成垃圾（边界全变 1、吸收失效）——复现本书的
    # 第一个坑。
    w = frac * N
    return np.exp(-(nsq ** n_order) / w ** (2 * n_order))


def ang_spec_multi_prop_vac(Uin, wvl, delta1, deltan, z, absorb='book'):
    """任意次部分传播的角谱法（书 Listing 8.1，式 8.18 的数值实现）。

    参数
        Uin    : U(r1)，源平面光场
        wvl    : 光波长 [m]
        delta1 : 源平面网格间距 δ1 [m]
        deltan : 观测平面网格间距 δn [m]
        z      : 中间/观测平面位置数组 [z2, z3, ..., zn]（不含 z1=0）[m]
        absorb : 'book'    —— 与 Listing 8.1 逐行一致（每段传播后都乘 sg，
                              观测面也被衰减，勘误 (c)）
                 'between' —— 只在中间平面衰减（式 8.9/8.16 的本意，严格版）
                 'none'    —— 不加吸收边界（算符链自洽性验证用）

    返回 (xn, yn, Uout)：观测面坐标与场。

    实现要点（与书代码逐行对应）：
        alpha = z/z(n)                  平面归一化位置（Table 8.2）
        delta = (1-alpha)*delta1 + alpha*deltan   网格间距线性插值（式 8.8）
        m = delta(2:n)./delta(1:n-1)    各段缩放因子（m_i = δi+1/δi）
        Q1 = Q[(1-m1)/Δz1, r1]          源面二次相位（先乘）
        循环 i = 1..n-1：
            Q2 = Q2[-Δzi/mi, fi]        频域二次相位
            U = sg .* ift2(Q2 .* ft2(U/mi, δi), δfi)   （式 8.18 的内层三算子）
        Q3 = Q[(m(n-1)-1)/(m(n-1)Δz(n-1)), rn]     观测面二次相位（后乘）
    """
    U = np.asarray(Uin, dtype=complex)
    N = U.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    sg = sg_absorbing_boundary(N)                 # 超高斯吸收边界（指数单位为样本）

    z = np.concatenate(([0.0], np.asarray(z, dtype=float)))
    nplanes = len(z)                              # 平面数 n
    Delta_z = np.diff(z)                          # Δzi = zi+1 - zi
    alpha = z / z[-1]                             # αi = zi/Δz
    delta = (1 - alpha) * delta1 + alpha * deltan  # 式 (8.8)：δi 线性插值
    m = delta[1:] / delta[:-1]                    # 缩放因子 mi

    x1, y1 = np.meshgrid(n * delta[0], n * delta[0])
    r1sq = x1 ** 2 + y1 ** 2
    Q1 = np.exp(1j * k / 2.0 * (1.0 - m[0]) / Delta_z[0] * r1sq)
    U = Q1 * U                                    # Q[(1-m1)/Δz1, r1]{U(r1)}
    # 注意：写成 Q1*U（而非 U*Q1）与第 6 章 Listing 6.5 的乘法顺序一致，
    # numpy 复数乘法的 FMA 收缩对操作数顺序不对称——顺序一致才能逐位复现。

    Z = Delta_z[-1]
    for idx in range(nplanes - 1):
        deltaf = 1.0 / (N * delta[idx])           # δfi = 1/(N δi)
        fX, fY = np.meshgrid(n * deltaf, n * deltaf)
        fsq = fX ** 2 + fY ** 2
        Z = Delta_z[idx]                          # 本段传播距离 Δzi
        Q2 = np.exp(-1j * np.pi ** 2 * 2.0 * Z / m[idx] / k * fsq)
        U = ift2(Q2 * ft2(U / m[idx], delta[idx]), deltaf)
        if absorb == 'book':
            U = sg * U                            # 勘误 (c)：含最后一段
        elif absorb == 'between' and idx < nplanes - 2:
            U = sg * U                            # 严格版：只衰减中间平面
        # absorb == 'none'：跳过

    xn, yn = np.meshgrid(n * delta[-1], n * delta[-1])
    rnsq = xn ** 2 + yn ** 2
    Q3 = np.exp(1j * k / 2.0 * (m[-1] - 1.0) / (m[-1] * Z) * rnsq)
    return xn, yn, Q3 * U


# ===== 式 (8.24)：部分传播的最大段距与最少段数 =====
def dz_max_partial(delta1, deltan, N, wvl):
    """Δzmax = min(δ1,δn)^2 · N / λ（式 8.24 右端，保守界限）。"""
    return min(delta1, deltan) ** 2 * N / wvl


def n_partial_book(Dz, dzmax):
    """书公式：至少 n = ceil(Δz/Δzmax) + 1 次部分传播（含 +1 裕量）。"""
    return int(np.ceil(Dz / dzmax)) + 1


def constraint4_N_exact_per_seg(delta1, deltan, wvl, Dz, nseg):
    """给定段数 nseg 时的精确逐段约束 4：max_i λΔzi/(δi δi+1)。

    返回 (Nmin_max, Delta_zi)：所有段中最苛刻的 N 下界与段距。
    δi 按式 (8.8) 线性插值（αi = zi/Δz，z 均匀分段时 zi = i·Δz/nseg）。
    """
    alpha = np.arange(nseg + 1) / nseg
    delta = (1 - alpha) * delta1 + alpha * deltan
    Dzi = Dz / nseg
    Nmin = wvl * Dzi / (delta[:-1] * delta[1:])
    return Nmin.max(), Dzi


def roi_metrics(xn, U, Uan, roi_half):
    """感兴趣区域 |x|<=roi_half、y=0 切片：max|ΔI|/峰值 与 相位残差。

    y=0 切片取中间行（MATLAB x2(N/2+1,:)；numpy 0-based 为第 N//2 行）。
    xn 的 x 坐标沿列变（meshgrid 'xy'），各行相同，取哪行都一样。
    """
    mid = U.shape[0] // 2
    sl = xn[mid, :]
    win = np.abs(sl) <= roi_half
    I, Ia = np.abs(U[mid, :]) ** 2, np.abs(Uan) ** 2
    pk = Ia.max()
    eI = np.max(np.abs(I[win] - Ia[win])) / pk
    mph = win & (Ia > 0.05 * pk)
    eph = np.max(np.abs(np.angle(U[mid, :][mph] / Uan[mph])))
    return eI, eph


def main(outdir=None):
    here = os.path.dirname(os.path.abspath(__file__))
    if outdir is None:
        outdir = here
    os.makedirs(outdir, exist_ok=True)
    results = []

    # ================= 验证 0：n=2 部分传播 == 单次角谱传播 =================
    print('===== 验证 0：n=2 部分传播（不加吸收边界）与 ang_spec_prop 逐位对比 =====')
    wvl = 1e-6
    N0 = 64
    d10, d20 = 20e-6, 55e-6
    Dz0 = 0.13
    rng = np.random.default_rng(8)
    U0 = rng.standard_normal((N0, N0)) + 1j * rng.standard_normal((N0, N0))
    _, _, U_multi = ang_spec_multi_prop_vac(U0, wvl, d10, d20, [Dz0], absorb='none')
    _, _, U_single = ang_spec_prop(U0, wvl, d10, d20, Dz0)
    dev0 = np.max(np.abs(U_multi - U_single)) / np.max(np.abs(U_single))
    print('    z=[%.2f m]（n=2 平面、1 段）：相对最大差 = %.2e' % (Dz0, dev0))
    assert dev0 == 0.0, 'n=2 部分传播应与单次角谱传播逐位一致'
    print('    通过：算符链 Q1/F/Q2/F-1/Q3 与第 6 章完全一致。')

    # ================= 验证 1：书例约束总账（式 8.24/8.25） =================
    print('===== 验证 1：书例约束总账（D1=2mm, λ=1µm, Δz=2m, 30 点跨孔径） =====')
    D1 = 2e-3
    d1 = D1 / 30                                    # δ1 = 66.67 µm
    variants = {'正文/图 8.5（D2=4 mm）': 4e-3, 'Listing 8.2（D2=6 mm）': 6e-3}
    pick = {}
    for name, D2 in variants.items():
        dn = D2 / 30
        c1 = constraint1_d2(d1, D1, D2, wvl, 2.0)
        c2N = constraint2_N(d1, dn, D1, D2, wvl, 2.0)
        c3lo, c3hi = constraint3_d2(d1, D1, wvl, 2.0, np.inf)
        zmax = dz_max_partial(d1, dn, 128, wvl)
        nbk = n_partial_book(2.0, zmax)
        pick[name] = (D2, dn, c1, c2N, c3hi, zmax, nbk)
        print('    [%s] δ1=%.1f µm, δn=%.1f µm' % (name, d1 * 1e6, dn * 1e6))
        print('        约束 1: δn <= %.0f µm；约束 2: N >= %.1f（N=128 %s）'
              % (c1 * 1e6, c2N, '满足' if 128 >= c2N else '轻微违反！'))
        print('        约束 3 (R=∞): δn <= %.0f µm；式 (8.25) Δzmax = %.4f m（书印 0.567 m）'
              % (c3hi * 1e6, zmax))
        print('        段数：书公式 ceil(2/Δzmax)+1 = %d 段；严格最小 ceil(2/Δzmax) = %d 段'
              % (nbk, int(np.ceil(2.0 / zmax))))
    results.append(('式(8.25) Δzmax [m]', zmax, 0.567, abs(zmax - 0.567)))

    # ================= 验证 2：Listing 8.2 复现（原书图 8.6） =================
    print('===== 验证 2：Listing 8.2 复现（N=128, 5 段 × 0.4 m） =====')
    z_planes = [0.4, 0.8, 1.2, 1.6, 2.0]            # 书： z = (1:5)*2/5
    out = {}
    for name, D2 in variants.items():
        dn = D2 / 30
        x1g, y1g = np.meshgrid((np.arange(128) - 64) * d1,
                               (np.arange(128) - 64) * d1)
        ap = rect(x1g / D1) * rect(y1g / D1)     # 与书一致：孔径边界半值
        xn, yn, Uout = ang_spec_multi_prop_vac(ap, wvl, d1, dn, z_planes)
        Dz = 2.0
        Uan = fresnel_prop_square_ap(xn[0, :], 0.0, D1, wvl, Dz)   # y=0 切片解析
        eI, eph = roi_metrics(xn, Uout, Uan, D2 / 2)
        out[name] = (xn, Uout, Uan, eI, eph)
        print('    [%s] δn=%.1f µm：ROI(±%.0f mm) |dI|/峰值 = %.2e，相位残差 = %.2e rad'
              % (name, dn * 1e6, D2 / 2 * 1e3, eI, eph))
    # 图 8.6 读数核对：采样点序列（判定 Listing 参数是实际运行版）
    xnL2 = out['Listing 8.2（D2=6 mm）'][0]
    UL = out['Listing 8.2（D2=6 mm）'][1]
    seq = [np.interp(k * 0.2e-3, xnL2[0, :],
                     np.abs(UL[UL.shape[0] // 2, :]) ** 2) for k in range(5)]
    print('    图 8.6 读数核对：I(k·0.2mm), k=0..4 = %s（原书图上 ≈ 2.55/2.25/1.55/0.8/0.35）'
          % '/'.join('%.2f' % v for v in seq))
    results.append(('L8.2 |dI|/pk', out['Listing 8.2（D2=6 mm）'][3], None, None))

    # ================= 验证 3：同一网格，Δz 加大——单次退化 vs 部分传播稳定 ==========
    print('===== 验证 3：固定网格（N=128, δ1=66.7µm, δn=200µm），Δz 加大 =====')
    dnL = 6e-3 / 30
    x1g, y1g = np.meshgrid((np.arange(128) - 64) * d1, (np.arange(128) - 64) * d1)
    ap = rect(x1g / D1) * rect(y1g / D1)         # 与书一致：孔径边界半值

    def roi_and_outer(xn, U, Uan, roi_half=3e-3):
        """ROI 误差 + 外环（卷绕敏感区）误差。"""
        mid = U.shape[0] // 2
        sl = xn[0, :]
        win = np.abs(sl) <= roi_half
        I, Ia = np.abs(U[mid, :]) ** 2, np.abs(Uan) ** 2
        pk = Ia.max()
        eI = np.max(np.abs(I[win] - Ia[win])) / pk
        mph = win & (Ia > 0.05 * pk)
        eph = np.max(np.abs(np.angle(U[mid, :][mph] / Uan[mph])))
        rout = (np.abs(sl) > roi_half) & (np.abs(sl) <= 0.8 * sl.max())
        eIout = np.max(np.abs(I[rout] - Ia[rout])) / pk
        return eI, eph, eIout

    print('    Dz[m]  单次约束4需N   单次: eI_ROI   eph     eI_外环 | 多段: 段数  eI_ROI   eph     eI_外环')
    scan = []
    for Dz in (2.0, 4.0, 6.0, 8.0):
        c4 = constraint4_N(d1, dnL, wvl, Dz)
        xnS, _, Us = ang_spec_prop(ap, wvl, d1, dnL, Dz)
        UaS = fresnel_prop_square_ap(xnS[0, :], 0.0, D1, wvl, Dz)
        e1, p1, o1 = roi_and_outer(xnS, Us, UaS)
        nseg = int(np.ceil(Dz / 0.4))            # 段距不超过书例的 0.4 m
        zpl = [Dz * i / nseg for i in range(1, nseg + 1)]
        xnM, _, Um = ang_spec_multi_prop_vac(ap, wvl, d1, dnL, zpl)
        UaM = fresnel_prop_square_ap(xnM[0, :], 0.0, D1, wvl, Dz)
        e2, p2, o2 = roi_and_outer(xnM, Um, UaM)
        scan.append((Dz, c4, e1, p1, o1, nseg, e2, p2, o2, xnS, Us, UaS, xnM, Um, UaM))
        print('    %4.0f   %6.0f        %.2e  %.2e  %.2e   |  %d      %.2e  %.2e  %.2e'
              % (Dz, c4, e1, p1, o1, nseg, e2, p2, o2))
    results.append(('Δz=8m 单次 |dI|/pk', scan[-1][2], None, None))
    results.append(('Δz=8m 20段 |dI|/pk', scan[-1][6], None, None))

    # ================= 验证 4：约束 3 与部分传播无关（式 8.20→8.22） =================
    print('===== 验证 4：约束 3 代换不变性（习题 2 的数值复核 + 有限 R 的补充） =====')
    print('    代数：式 (8.21) 代入 δ2=(1-α2)δ1+α2δn、Δz1=α2Δz，消去公因子 α2 后')
    print('    与式 (8.22) 恒等（R=∞ 时严格；推导见笔记习题 2 拓展）。')
    worst = 0.0
    for _ in range(200):
        d1r = 10 ** rng.uniform(-5, -3)
        dnr = 10 ** rng.uniform(-5, -3)
        Dzr = 10 ** rng.uniform(-1, 1)
        wvlr = 10 ** rng.uniform(-7, -6)
        Rr = 10 ** rng.uniform(-1, 3) * np.sign(rng.uniform(-1, 1) + 0.5)
        a2 = rng.uniform(0.05, 0.95)
        d2r = (1 - a2) * d1r + a2 * dnr          # 式 (8.8)：δ2
        lo21 = (1 + a2 * Dzr / Rr) * d1r - wvlr * a2 * Dzr / d1r   # 式 (8.21) 左（对 δ2）
        hi21 = (1 + a2 * Dzr / Rr) * d1r + wvlr * a2 * Dzr / d1r   # 式 (8.21) 右
        lo22 = (1 + Dzr / Rr) * d1r - wvlr * Dzr / d1r             # 式 (8.22) 左（对 δn）
        hi22 = (1 + Dzr / Rr) * d1r + wvlr * Dzr / d1r             # 式 (8.22) 右
        # 代数等价：δ2 <= hi21 ⇔ δn <= hi22（两边消去公因子 α2）。浮点核算：
        # 界可能过零，分母用区间半宽（2λΔz/D1，恒正）。
        halfw22 = hi22 - lo22
        worst = max(worst, abs((lo21 - (1 - a2) * d1r) / a2 - lo22) / halfw22,
                    abs((hi21 - (1 - a2) * d1r) / a2 - hi22) / halfw22)
        eps = 1e-9                                # 抗刀刃 cases 的相对裕量
        ok21 = (lo21 - eps * abs(lo21)) <= d2r <= (hi21 + eps * abs(hi21))
        ok22 = (lo22 - eps * abs(lo22)) <= dnr <= (hi22 + eps * abs(hi22))
        assert ok21 == ok22
    print('    200 组随机参数：式 (8.21) 界折算到 δn 空间与式 (8.22) 界的最大差'
          '（以区间半宽为单位）= %.2e；"δ2 是否满足 (8.21)" 与 "δn 是否满足'
          ' (8.22)" 判定 200/200 一致' % worst)
    print('    有限 R 补充（本笔记推导）：第 i 段（i>1）的约束 3 中心从')
    print('    (1+Δz/R)δ1 偏移到 δ1+Δz·δi/R（用第 i 平面的束径 D_i·曲率 R+z_i 重推），')
    print('    R→∞ 时偏移归零——书例 R=∞ 不受影响；正文"与段数无关"对 R=∞ 严格成立。')

    # ================= 验证 5：吸收边界施加位置的勘误 (c) 实锤 =================
    print('===== 验证 5：sg 施加位置——book 模式 vs between 模式 =====')
    _, _, U_book = ang_spec_multi_prop_vac(ap, wvl, d1, dnL, z_planes, absorb='book')
    _, _, U_btwn = ang_spec_multi_prop_vac(ap, wvl, d1, dnL, z_planes, absorb='between')
    Uan5 = out['Listing 8.2（D2=6 mm）'][2]
    edge_r = 60                                     # 距中心 60 个样本（边缘区）
    amp_ratio = np.abs(U_book[64, 64 + edge_r]) / np.abs(U_btwn[64, 64 + edge_r])
    print('    ROI 内两种模式 |dI|/峰值 相对解析分别 %.2e / %.2e（sg≈1，无差别）'
          % (roi_metrics(xnL2, U_book, Uan5, 3e-3)[0],
             roi_metrics(xnL2, U_btwn, Uan5, 3e-3)[0]))
    print('    网格边缘（r=60 样本）振幅比 |U_book|/|U_between| = %.3f'
          % amp_ratio)
    print('    —— book 模式在观测面多衰减一次（勘误 c 实锤）；ROI 内 sg≈1 两者一致。')

    # ================= 验证 6：吸收边界函数实现核对（式 8.1 / Listing 8.1） =================
    print('===== 验证 6：超高斯吸收边界实现核对 =====')
    sg = sg_absorbing_boundary(128)
    nvec = (np.arange(128) - 64).astype(float)      # float！整数 8 次幂会 int64 溢出
    NX, NY = np.meshgrid(nvec, nvec)
    nsq = NX ** 2 + NY ** 2
    w = 0.47 * 128
    sg_ref = np.exp(-nsq ** 8 / w ** 16)            # 书 Listing 8.1 原式
    dev6 = np.max(np.abs(sg - sg_ref))
    # 式 (8.1)：g_sg = exp(-(r/σ)^n)，σ = 0.47N（样本）、n = 16
    sg_e81 = np.exp(-(np.sqrt(nsq) / (0.47 * 128)) ** 16)
    dev81 = np.max(np.abs(sg - sg_e81))
    print('    与 Listing 8.1 原式 exp(-nsq.^8/w^16) 最大差 = %.2e' % dev6)
    print('    与式 (8.1) 形式 exp(-(r/σ)^n)（σ=0.47N, n=16）最大差 = %.2e' % dev81)
    print('    边界值抽查：中心 sg=%.3f；r=0.47N 处 %.3f；网格角 r=%.1f 处 %.2e'
          % (sg[64, 64], sg[64, 64 + 60],
             np.sqrt(2) * 64, sg[0, 0]))
    assert dev6 == 0.0 and dev81 < 1e-14            # (r/σ)^16 与 r^16/σ^16 差 1 ULP

    # ================= 图 8-7：书例复现（约束图 + 辐照度/相位） =================
    fig7, axs7 = plt.subplots(1, 3, figsize=(15.5, 4.8))
    ax = axs7[0]
    d1v = np.linspace(2e-6, 260e-6, 400)
    dnv = np.linspace(5e-6, 550e-6, 400)
    D1g, Dng = np.meshgrid(d1v, dnv)
    log2N = np.log2(np.maximum(constraint2_N(D1g, Dng, 4e-3, 4e-3, wvl, 2.0), 1e-3))
    cs = ax.contourf(D1g * 1e6, Dng * 1e6, log2N,
                     levels=np.arange(3, 14.5, 0.5), cmap='Greys')
    clb = fig7.colorbar(cs, ax=ax)
    clb.set_label(r'$\log_2 N$（约束 2，式 7.20）')
    c1v = constraint1_d2(d1v, 2e-3, 4e-3, wvl, 2.0)          # 正文参数（截距 1000 µm）
    ax.plot(d1v * 1e6, np.clip(c1v, 0, 550) * 1e6, 'k--', lw=2,
            label='约束 1（式 7.14，截距 1000 $\mu$m）')
    ax.plot(d1v * 1e6, np.clip(500e-6 - 2 * d1v, 0, 550) * 1e6, color='0.45',
            ls=':', lw=2, label='图 8.5 实测虚线（截距 500 $\mu$m，勘误 b）')
    ax.plot(d1 * 1e6, 200e-6 * 1e6, 'g*', ms=14, zorder=5)
    ax.annotate('Listing 8.2\n(66.7, 200)', (d1 * 1e6, 200e-6 * 1e6),
                textcoords='offset points',
                xytext=(8, 6), fontsize=8, color='g')
    ax.plot(d1 * 1e6, 4e-3 / 30 * 1e6, 'o', mfc='none', mec='C1', ms=10, mew=2, zorder=5)
    ax.annotate('正文/图 8.5 白叉\n(66.7, 133)', (d1 * 1e6, 133.3),
                textcoords='offset points', xytext=(10, -14), fontsize=8, color='C1')
    ax.set_xlabel(r'$\delta_1$ [$\mu$m]'); ax.set_ylabel(r'$\delta_n$ [$\mu$m]')
    ax.set_xlim(0, 260); ax.set_ylim(0, 550)
    ax.set_title('(a) 书例约束图（等高线按 D2=4mm、Δz=2m 重建）', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    xnL_, UL, UanL, eIL, ephL = out['Listing 8.2（D2=6 mm）']
    mid = UL.shape[0] // 2
    sl = xnL_[0, :] * 1e3
    ax = axs7[1]
    ax.plot(sl, np.abs(UanL) ** 2 / np.abs(UanL).max(), 'k-', lw=1.2, label='解析')
    ax.plot(sl, np.abs(UL[mid, :]) ** 2 / np.abs(UanL).max(), 'rx', ms=4,
            mew=1.2, label='数值（N=128, 5 段）')
    ax.set_xlim(-3, 3)
    ax.set_xlabel(r'$x_n$ [mm]'); ax.set_ylabel('归一化辐照度')
    ax.set_title(r'(b) Listing 8.2 复现：ROI 辐照度（|dI|/pk = %.1e）' % eIL, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs7[2]
    mph = (np.abs(xnL_[0, :]) <= 3e-3) & (np.abs(UanL) > 0.05 * np.abs(UanL).max())
    ax.plot(sl, np.angle(UanL), 'k-', lw=1.2, label='解析')
    ax.plot(sl[mph], np.angle(UL[mid, :][mph]), 'rx', ms=4, mew=1.2, label='数值')
    ax.set_xlim(-3, 3)
    ax.set_xlabel(r'$x_n$ [mm]'); ax.set_ylabel('相位 [rad]')
    ax.set_title(r'(c) ROI 相位（残差 %.1e rad，|U|>5%% 峰值）' % ephL, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig7.suptitle('图 8-7（自绘仿真）工程 1：Listing 8.2 复现（D1=2 mm, λ=1 $\mu$m, Δz=2 m, N=128, 5 段×0.4 m）'
                  '（(a) 中约束 1 按正文参数画：截距 1000 $\mu$m；原书图 8.5 虚线截距实测 500 $\mu$m）',
                  fontsize=12)
    fig7.tight_layout(rect=[0, 0, 1, 0.9])
    fp7 = os.path.join(outdir, 'fig-8-7.png')
    fig7.savefig(fp7, dpi=150)
    print('    已保存 %s' % fp7)

    # ================= 图 8-8：固定网格加大 Δz——单次退化 vs 部分传播稳定 ==========
    fig8, axs8 = plt.subplots(1, 3, figsize=(15.5, 4.8))
    s8 = scan[-1]                                 # Δz = 8 m
    Dz8, c48, xnS8, Us8, UaS8, xnM8, Um8, UaM8 = (s8[0], s8[1], s8[9], s8[10],
                                                  s8[11], s8[12], s8[13], s8[14])
    mid8 = Us8.shape[0] // 2
    ax = axs8[0]
    ax.plot(xnS8[0, :] * 1e3, np.abs(UaS8) ** 2 / np.abs(UaS8).max(), 'k-', lw=1.2,
            label='解析')
    ax.plot(xnS8[0, :] * 1e3, np.abs(Us8[mid8, :]) ** 2 / np.abs(UaS8).max(),
            'C3+-', ms=4, lw=0.8,
            label=r'单次传播（N=128，约束 4 需 N$\geq$%d）' % c48)
    ax.plot(xnM8[0, :] * 1e3, np.abs(Um8[mid8, :]) ** 2 / np.abs(UaM8).max(),
            'b.', ms=4, label='20 段×0.4 m 部分传播（N=128）')
    ax.set_xlim(-3, 3)
    ax.set_xlabel(r'$x_n$ [mm]'); ax.set_ylabel('归一化辐照度')
    ax.set_title(r'(a) Δz=8 m、ROI 辐照度：单次误差 %.1e vs 多段 %.1e'
                 % (s8[2], s8[6]), fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ax = axs8[1]
    mph8 = (np.abs(xnS8[0, :]) <= 3e-3) & (np.abs(UaS8) > 0.05 * np.abs(UaS8).max())
    ax.plot(xnS8[0, :] * 1e3, np.angle(UaS8), 'k-', lw=1.2, label='解析')
    ax.plot(xnS8[0, :][mph8] * 1e3, np.angle(Us8[mid8, :][mph8]), 'C3+', ms=5, mew=1.2,
            label='单次传播')
    ax.plot(xnM8[0, :][mph8] * 1e3, np.angle(Um8[mid8, :][mph8]), 'b.', ms=4,
            label='20 段部分传播')
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3.3, 3.3)
    ax.set_xlabel(r'$x_n$ [mm]'); ax.set_ylabel('相位 [rad]')
    ax.set_title('(b) ROI 相位：单次残差 %.1e rad / 多段 %.1e rad' % (s8[3], s8[7]),
                 fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ax = axs8[2]
    dzs = [s[0] for s in scan]
    ax.semilogy(dzs, [s[2] for s in scan], 'C3o-', lw=1.4, label='单次传播：ROI |dI|/pk')
    ax.semilogy(dzs, [s[6] for s in scan], 'bo-', lw=1.4, mfc='none',
                label='多段部分传播：ROI |dI|/pk')
    ax.semilogy(dzs, [s[3] for s in scan], 'C3^--', lw=1.0, ms=5, alpha=0.6,
                label='单次：ROI 相位残差')
    ax.semilogy(dzs, [s[7] for s in scan], 'b^--', lw=1.0, ms=5, alpha=0.6,
                label='多段：ROI 相位残差')
    for s in scan:
        ax.annotate('N$\\geq$%d' % s[1], (s[0], s[2]), textcoords='offset points',
                    xytext=(2, 7), fontsize=7, color='C3')
    ax.set_xlabel(r'总传播距离 $\Delta z$ [m]')
    ax.set_ylabel('误差（ROI，y=0 切片）')
    ax.set_title('(c) 固定网格加大距离：单次退化、多段稳定\n'
                 r'（段数 = ceil($\Delta z$/0.4 m)，约束 4 逐段满足）', fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3, which='both')
    fig8.suptitle('图 8-8（自绘仿真）工程 1：约束 4 是部分传播要消灭的敌人'
                  '（D1=2 mm, δ1=66.7 $\mu$m, δn=200 $\mu$m, N=128, λ=1 $\mu$m；'
                  '段距 ≤0.4 m 时每段约束 4 只需 N≤75）', fontsize=12)
    fig8.tight_layout(rect=[0, 0, 1, 0.9])
    fp8 = os.path.join(outdir, 'fig-8-8.png')
    fig8.savefig(fp8, dpi=150)
    print('    已保存 %s' % fp8)

    # ================= 结果总账 =================
    print('===== 结果总账 =====')
    for name, val, ref, dev in results:
        if ref is None:
            print('    %s = %s' % (name, val))
        else:
            print('    %s = %s（书印 %s，差 %s）' % (name, val, ref, dev))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
