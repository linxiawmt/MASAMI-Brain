# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 6 章（工程 1）=====
真空中的菲涅耳衍射：单步 / 两步 / 角谱 / ABCD 四种传播器
（原书 6.1-6.5 节，Listings 6.1-6.7）

原书要点：
    Fresnel 积分（式 6.1）的两种数值形态：
      单步 FT 形态（式 6.5，Listing 6.1）：U = Q V F Q —— 一次 FT，
        观测面采样 delta2 = lambda*Dz/(N*delta1) 固定（式 6.16）；
      卷积形态（式 6.6）：两步传播（式 6.17-6.30，Listing 6.3）与
        角谱传播（式 6.31/6.50，Listing 6.5）都引入缩放参数 m = delta2/delta1。
    算子记号（式 6.7-6.13）：Q/V/F/F^-1/R/Q2，贯穿第 6-8 章。
    简单光学系统（式 6.68-6.81）：ABCD 阵 + 广义 Huygens-Fresnel 积分
    （Collins，式 6.77）的 Lambert-Fraser 卷积化（式 6.79-6.81）；
    Listing 6.6 实现其自由 m 版本（本笔记补出算子链，见验证 5）。

数值验证（运行本脚本全部自动打印）：
    (0) ft2/ift2 往返误差 + Parseval；
    (1) 书例复现（Listing 6.2，N=1024 方形孔径 D=2mm, Dz=1m）：one_step vs
        Fresnel 积分解析闭式，辐照度/相位切片对照（原书图 6.1）；
    (2) 书例复现（Listing 6.4，m=10.24）：two_step vs one_step vs 解析
        （原书图 6.4）；表 6.2 数值核对；式 (6.19)-(6.23)、(6.55)/(6.56) 核对；
    (3) 非对称场（离轴矩形+圆盘，m=2.5）：
        - 书印 Listing 6.3 "-" 支的镜像效应（对称测试场不可见！）：
          输出 = 真值的镜像（|U| 逐位一致），翻转中间场修正后消失；
        - 修正 "-" 支 two_step_fixed 与 ang_spec 逐位一致（1e-14，
          全部测试的 m 成立）——两条链在离散层面是同一算子；
        - "+" 支 / 式(6.66)修正链是另一条等价离散路径（差 3~5e-2，
          随 m->1 与网格加密收敛，细网格仲裁实测）；
    (4) 式 (6.66) 勘误验证：书印链 vs 修正链（= 式 6.67 顶符号）——
        实输入场下书印链输出 = -conj(真值)（逐点、机器精度 4e-16）；
        点源解析检验：书印链给 e^{-i k|r2-rc|^2/(2Dz)}（相位反号，见笔记）；
    (5) Listing 6.7 复现：ABCD（f=inf 自由空间）== ang_spec（机器精度）；
        ABCD（薄透镜 f=2m + 自由空间）== one_step(孔径 x 发散相位) ——
        Collins 积分与"薄透镜相位 + Fresnel 传播"的等价性核对；
    (6) 采样判据核对：源面/观测面二次相位的局部空间频率 vs Nyquist。

原书勘误（本脚本按修正实现并数值验证，详见笔记）：
    (a) 式 (6.66) 的算子链不是 (6.50) 的正确 m->-m 代换：Q 系数应为
        (m+1)/(m*Dz) 与 (1+m)/Dz（书印 (1-m)/(m*Dz) 与 -(1-m)/Dz），
        F^-1 的求值点应为 -r2/m（书印 +r2/m）；(6.67) 顶符号才是对的。
    (b) Listing 6.3 对负距离步用 |Dz1| 建中间面坐标网格，坐标对应翻转：
        m != 1 时 "-" 支输出是真值的镜像（对称测试场不可见）。

MATLAB/Python 差异（对应 .m 文件内注释）：
    - numpy 0-based / MATLAB 1-based：居中排列的零坐标样本位置不同；
    - np.meshgrid 默认 'xy' 与 MATLAB meshgrid 一致；
    - 偶数网格的居中坐标反转 P(-r)：两语言都用 flip + roll 1
      （MATLAB: circshift(rot90(P,2),[1 1])；numpy: roll(flip(P),(1,1))）；
    - scipy.special.fresnel 直接给 Fresnel 积分；MATLAB 无内置（符号工具箱
      除外），.m 里用 Abramowitz-Stegun 7.3.32/33 有理近似（误差 < 3e-7）。
用法：python ch6_fresnel_prop.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib + scipy；base.py（rect/circle/ft2_centered，同仓库
      python/fourier/common）。
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import fresnel

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ---- 复用《计算傅里叶光学》笔记的公共基础模块（同仓库） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import rect, circle, ft2_centered as ft2  # noqa: E402


# ===== 书 Listing 2.6：ift2（定标正确的二维逆 DFT） =====
def ift2(G, delta_f):
    """连续逆傅里叶变换的 DFT 近似（居中输出）：
        g(r) ~= (N*delta_f)^2 * IFFT{G}
    numpy.fft.ifft2 自带 1/N^2，乘 (N*delta_f)^2 回到连续量纲。"""
    N = G.shape[0]
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(G))) * (N * delta_f) ** 2


def rollflip(U):
    """偶数网格居中坐标反转：(Mf)(r) = f(-r) 的居中采样。
    numpy: roll(flip(U), (1,1))；MATLAB: circshift(rot90(U,2), [1 1])。"""
    return np.roll(np.rot90(U, 2), (1, 1), axis=(0, 1))


# ===== 书 Listing 6.1：one_step_prop =====
def one_step_prop(Uin, wvl, d1, Dz):
    """单步 Fresnel 传播（式 6.5；全局相位 e^{ik*Dz} 按原书代码约定省略）。

    返回 (x2, y2, Uout)。观测面采样 delta2 = lambda*Dz/(N*d1)（式 6.16），
    无法由用户指定 —— 这正是两步/角谱法的动机。"""
    N = Uin.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    d2 = wvl * Dz / (N * d1)                       # 式 (6.16)
    x2, y2 = np.meshgrid(n * d2, n * d2)
    Uout = (1.0 / (1j * wvl * Dz)
            * np.exp(1j * k / (2.0 * Dz) * (x2 ** 2 + y2 ** 2))
            * ft2(Uin * np.exp(1j * k / (2.0 * Dz) * (x1 ** 2 + y1 ** 2)), d1))
    return x2, y2, Uout


# ===== 书 Listing 6.3：two_step_prop（书印 "-" 支） =====
def two_step_prop(Uin, wvl, d1, d2, Dz):
    """两步 Fresnel 传播，书印 Listing 6.3 的 "-" 支：Dz1 = Dz/(1-m)（式 6.25）。

    注意（勘误 b）：Dz1（m>1）或 Dz2（m<1）必有一个为负；原书对负距离步用
    |Dz| 建坐标网格（第 14 行 d1a = wvl*abs(Dz1)/(N*d1)），坐标对应翻转一次，
    输出为真值的镜像（对称测试场不可见）。本函数忠实复现书印行为；
    修正版见 two_step_prop_fixed / two_step_prop_plus。"""
    N = Uin.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    m = d2 / d1
    Dz1 = Dz / (1.0 - m)
    d1a = wvl * abs(Dz1) / (N * d1)                # 式 (6.19)（书印取 |Dz1|）
    x1a, y1a = np.meshgrid(n * d1a, n * d1a)
    Uitm = (1.0 / (1j * wvl * Dz1)
            * np.exp(1j * k / (2.0 * Dz1) * (x1a ** 2 + y1a ** 2))
            * ft2(Uin * np.exp(1j * k / (2.0 * Dz1) * (x1 ** 2 + y1 ** 2)), d1))
    Dz2 = Dz - Dz1                                 # Dz2 = -m*Dz1
    x2, y2 = np.meshgrid(n * d2, n * d2)
    Uout = (1.0 / (1j * wvl * Dz2)
            * np.exp(1j * k / (2.0 * Dz2) * (x2 ** 2 + y2 ** 2))
            * ft2(Uitm * np.exp(1j * k / (2.0 * Dz2) * (x1a ** 2 + y1a ** 2)),
                  d1a))
    return x2, y2, Uout


def two_step_prop_fixed(Uin, wvl, d1, d2, Dz):
    """两步传播 "-" 支的修正版：负距离步的输出按居中反转重排，
    消除书印实现的镜像（其余与书印逐行对应）。
    数值事实：修正后与 ang_spec 逐位一致（1e-14，全部测试的 m）。"""
    N = Uin.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    m = d2 / d1
    Dz1 = Dz / (1.0 - m)
    d1a = wvl * abs(Dz1) / (N * d1)
    x1a, y1a = np.meshgrid(n * d1a, n * d1a)
    Uitm = (1.0 / (1j * wvl * Dz1)
            * np.exp(1j * k / (2.0 * Dz1) * (x1a ** 2 + y1a ** 2))
            * ft2(Uin * np.exp(1j * k / (2.0 * Dz1) * (x1 ** 2 + y1 ** 2)), d1))
    if Dz1 < 0:                                    # 修正：负距离步输出镜像翻转
        Uitm = rollflip(Uitm)
    Dz2 = Dz - Dz1
    x2, y2 = np.meshgrid(n * d2, n * d2)
    Uout = (1.0 / (1j * wvl * Dz2)
            * np.exp(1j * k / (2.0 * Dz2) * (x2 ** 2 + y2 ** 2))
            * ft2(Uitm * np.exp(1j * k / (2.0 * Dz2) * (x1a ** 2 + y1a ** 2)),
                  d1a))
    if Dz2 < 0:                                    # m<1 时第二步为负距离
        Uout = rollflip(Uout)
    return x2, y2, Uout


def two_step_prop_plus(Uin, wvl, d1, d2, Dz):
    """两步传播 "+" 支（式 6.25 的 "+"，图 6.3 几何）：中间面夹在源面与
    观测面之间，Dz1 = Dz/(1+m) > 0、Dz2 = m*Dz/(1+m) > 0（m>0），
    两个步都是正距离，无镜像问题。与 "-" 支在连续域严格等价
    （|Dz2/Dz1| = m，式 6.24/6.30）；离散路径不同（见验证 3 的收敛检查）。"""
    N = Uin.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    m = d2 / d1
    Dz1 = Dz / (1.0 + m)
    d1a = wvl * abs(Dz1) / (N * d1)
    x1a, y1a = np.meshgrid(n * d1a, n * d1a)
    Uitm = (1.0 / (1j * wvl * Dz1)
            * np.exp(1j * k / (2.0 * Dz1) * (x1a ** 2 + y1a ** 2))
            * ft2(Uin * np.exp(1j * k / (2.0 * Dz1) * (x1 ** 2 + y1 ** 2)), d1))
    Dz2 = Dz - Dz1
    x2, y2 = np.meshgrid(n * d2, n * d2)
    Uout = (1.0 / (1j * wvl * Dz2)
            * np.exp(1j * k / (2.0 * Dz2) * (x2 ** 2 + y2 ** 2))
            * ft2(Uitm * np.exp(1j * k / (2.0 * Dz2) * (x1a ** 2 + y1a ** 2)),
                  d1a))
    return x2, y2, Uout


# ===== 书 Listing 6.5：ang_spec_prop =====
def ang_spec_prop(Uin, wvl, d1, d2, Dz):
    """角谱传播（式 6.50 的 m 分支，Listing 6.5）：
        Uout = Q3 * ift2( Q2 * ft2( Q1 * Uin / m ) )
    Q1 = Q[(1-m)/Dz]（空间网格）、Q2 = Q2[-Dz/m]（频率网格，
    = e^{-i*pi*lambda*Dz*f^2/m}）、Q3 = Q[(m-1)/(m*Dz)]（空间网格）。
    delta2 = m*delta1 由用户自由定标。"""
    N = Uin.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    r1sq = x1 ** 2 + y1 ** 2
    df1 = 1.0 / (N * d1)
    fX, fY = np.meshgrid(n * df1, n * df1)
    fsq = fX ** 2 + fY ** 2
    m = d2 / d1
    x2, y2 = np.meshgrid(n * d2, n * d2)
    r2sq = x2 ** 2 + y2 ** 2
    Q1 = np.exp(1j * k / 2.0 * (1.0 - m) / Dz * r1sq)
    Q2 = np.exp(-1j * np.pi ** 2 * 2.0 * Dz / m / k * fsq)
    Q3 = np.exp(1j * k / 2.0 * (m - 1.0) / (m * Dz) * r2sq)
    Uout = Q3 * ift2(Q2 * ft2(Q1 * Uin / m, d1), df1)
    return x2, y2, Uout


# ===== 验证 4：式 (6.66) 的两种实现（勘误对照） =====
def ang_spec_mbranch(Uin, wvl, d1, d2, Dz, corrected=True):
    """式 (6.66)（m' = -m 分支）的算子链实现。

    corrected=True（按 m'=-m 正确代换，即式 (6.67) 顶符号）：
        U = Q[(m+1)/(m*Dz)] * F^-1[f1, -r2/m] * Q2[+Dz/m] * F * Q[(1+m)/Dz] * (-1/m)
    corrected=False（按书印 (6.66)）：
        U = Q[(1-m)/(m*Dz)] * F^-1[f1, +r2/m] * Q2[+Dz/m] * F * Q[-(1-m)/Dz] * (-1/m)
    实输入场下书印链输出 = -conj(真值)（逐点解析恒等，验证 4 实测 4e-16）。
    """
    N = Uin.shape[0]
    k = 2.0 * np.pi / wvl
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    r1sq = x1 ** 2 + y1 ** 2
    df1 = 1.0 / (N * d1)
    fX, fY = np.meshgrid(n * df1, n * df1)
    fsq = fX ** 2 + fY ** 2
    m = d2 / d1
    x2, y2 = np.meshgrid(n * d2, n * d2)
    r2sq = x2 ** 2 + y2 ** 2
    if corrected:
        Q1 = np.exp(1j * k / 2.0 * (1.0 + m) / Dz * r1sq)
        Q3 = np.exp(1j * k / 2.0 * (m + 1.0) / (m * Dz) * r2sq)
    else:
        Q1 = np.exp(-1j * k / 2.0 * (1.0 - m) / Dz * r1sq)       # 书印 Q[-(1-m)/Dz]
        Q3 = np.exp(1j * k / 2.0 * (1.0 - m) / (m * Dz) * r2sq)  # 书印 Q[-(m-1)/(m*Dz)]
    Q2 = np.exp(1j * np.pi ** 2 * 2.0 * Dz / m / k * fsq)        # Q2[+Dz/m]
    W = (-1.0 / m) * ift2(Q2 * ft2(Q1 * Uin, d1), df1)           # 本征 IFT 网格
    if corrected:
        W = rollflip(W)          # F^-1[f1, -r2/m]：本征网格居中反转采样
    return x2, y2, Q3 * W


# ===== 书 Listing 6.6：ang_spec_propABCD =====
def ang_spec_propABCD(Uin, wvl, d1, d2, ABCD):
    """角谱 + ABCD 光线矩阵（式 6.79-6.81 的自由 m 版本，Listing 6.6）。

    算子链（本笔记补出的推导；自由空间 A=D=1 时退化为式 6.50）：
        Uout = Q[(D-1/m)/B] * ift2( Q2[-B/m] * ft2( Q[(A-m)/B] * Uin / m ) )
    C 不显式出现（AD-BC=1 保证 (6.79) 自洽）。
    原书签名注释把 wvl 印成 wwl（笔误），此处按实现更正。"""
    N = Uin.shape[0]
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    r1sq = x1 ** 2 + y1 ** 2
    df1 = 1.0 / (N * d1)
    fX, fY = np.meshgrid(n * df1, n * df1)
    fsq = fX ** 2 + fY ** 2
    m = d2 / d1
    x2, y2 = np.meshgrid(n * d2, n * d2)
    r2sq = x2 ** 2 + y2 ** 2
    A = ABCD[0, 0]
    B = ABCD[0, 1]
    D = ABCD[1, 1]
    Q1 = np.exp(1j * np.pi / (wvl * B) * (A - m) * r1sq)
    Q2 = np.exp(-1j * np.pi * wvl * B / m * fsq)
    Q3 = np.exp(1j * np.pi / (wvl * B) * (D - 1.0 / m) * r2sq)
    Uout = Q3 * ift2(Q2 * ft2(Q1 * Uin / m, d1), df1)
    return x2, y2, Uout


# ===== 方形孔径 Fresnel 衍射解析解（原书 fresnel_prop_square_ap） =====
def fresnel_prop_square_ap(x2, y2, D, wvl, Dz):
    """边长 D 的单位振幅方形孔径的 Fresnel 衍射解析解（式 6.1 可分离积分）。

        U(x2,y2) = (e^{ik*Dz}/(i*lambda*Dz)) * I(x2) * I(y2)
        I(x) = int_{x-D/2}^{x+D/2} exp(i*pi*u^2/(lambda*Dz)) du
             = sqrt(lambda*Dz/2) * {[C(t2)-C(t1)] + i[S(t2)-S(t1)]}
        t = sqrt(2/(lambda*Dz)) * u   （scipy.special.fresnel 的约定）
    全局相位约定：与原书代码一致，省略 e^{ik*Dz}（见笔记 6.1 节注记）。
    """
    a = 1.0 / (wvl * Dz)                 # 相位 = pi*a*u^2
    s = np.sqrt(2.0 * a)

    def I(x):
        S1, C1 = fresnel(s * (x - D / 2.0))
        S2, C2 = fresnel(s * (x + D / 2.0))
        return np.sqrt(1.0 / (2.0 * a)) * ((C2 - C1) + 1j * (S2 - S1))

    return I(x2) * I(y2) / (1j * wvl * Dz)


def main(outdir=None):
    # ---- 验证 0：ft2/ift2 往返 + Parseval ----
    N0, L0 = 64, 1.0
    d0 = L0 / N0
    rng = np.random.default_rng(0)
    g0 = rng.standard_normal((N0, N0)) + 1j * rng.standard_normal((N0, N0))
    df0 = 1.0 / (N0 * d0)
    rt = ift2(ft2(g0, d0), df0)
    Eg = np.sum(np.abs(g0) ** 2) * d0 ** 2
    EG = np.sum(np.abs(ft2(g0, d0)) ** 2) * df0 ** 2
    print('===== 验证 0：ft2/ift2 往返 + Parseval =====')
    print('    max|ift2(ft2(g)) - g| = %.2e' % np.max(np.abs(rt - g0)))
    print('    Parseval: sum|g|^2*d^2 = %.6e vs sum|G|^2*df^2 = %.6e（相对差 %.1e）'
          % (Eg, EG, abs(Eg - EG) / Eg))

    # ---- 原书 Listing 6.2/6.4 的参数（方形孔径，Dz = 1 m，N = 1024） ----
    N = 1024                     # 原书 Listing 6.2 的规模
    L = 1e-2
    d1 = L / N
    D = 2e-3
    wvl = 1e-6
    k = 2.0 * np.pi / wvl
    Dz = 1.0
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    ap = rect(x1 / D) * rect(y1 / D)

    print('===== 验证 1：Listing 6.2 单步传播 vs 解析（原书图 6.1） =====')
    x2, y2, U1 = one_step_prop(ap, wvl, d1, Dz)
    Uan = fresnel_prop_square_ap(x2, y2, D, wvl, Dz)
    win = (np.abs(x2) <= 3e-3) & (np.abs(y2) <= 3e-3)
    I1, Ian = np.abs(U1) ** 2, np.abs(Uan) ** 2
    pkA = np.max(Ian[win])
    strong = win & (Ian > 0.05 * pkA)              # 避开解析零点（比值发散）
    ratio = U1[strong] / Uan[strong]
    ph_res = np.abs(np.angle(ratio) - np.mean(np.angle(ratio)))
    print('    delta2 = %.4f mm（式 6.16），观测网格半宽 %.1f mm'
          % (wvl * Dz / (N * d1) * 1e3, np.max(x2) * 1e3))
    print('    轴上峰值辐照度：数值 %.4f vs 解析 %.4f（相对差 %.1e）'
          % (I1[N // 2, N // 2], Ian[N // 2, N // 2],
             abs(I1[N // 2, N // 2] - Ian[N // 2, N // 2]) / Ian[N // 2, N // 2]))
    print('    ±3 mm 窗内：max|U1-Uan|/max|Uan| = %.2e；'
          '辐照度最大绝对差/峰值 = %.2e'
          % (np.max(np.abs(U1[win] - Uan[win])) / np.max(np.abs(Uan[win])),
             np.max(np.abs(I1[win] - Ian[win])) / pkA))
    print('    |U|>5%% 峰值处：复场比值 |比值-1| max = %.2e；'
          '相位残差(去均值) max = %.2e rad（离散化水平）'
          % (np.max(np.abs(np.abs(ratio) - 1.0)), np.max(ph_res)))

    print('===== 验证 2：Listing 6.4 两步传播（m=10.24）+ 表 6.2 + 网格关系 =====')
    d2_book = wvl * Dz / (N * d1)          # 原书取 one-step 的 delta2 -> m = 10.24
    m_book = d2_book / d1
    x2b, y2b, U2 = two_step_prop(ap, wvl, d1, d2_book, Dz)
    x2p, y2p, U2p = two_step_prop_plus(ap, wvl, d1, d2_book, Dz)
    I2 = np.abs(U2) ** 2
    print('    m = delta2/delta1 = %.2f：Dz1 = Dz/(1-m) = %.4f m（<0，图 6.2 几何）'
          % (m_book, Dz / (1 - m_book)))
    print('    轴上峰值辐照度：two_step %.4f vs 解析 %.4f（相对差 %.1e）'
          % (I2[N // 2, N // 2], Ian[N // 2, N // 2],
             abs(I2[N // 2, N // 2] - Ian[N // 2, N // 2]) / Ian[N // 2, N // 2]))
    print('    two_step(书印"-"支) vs one_step 复场 max|diff|/max|U| = %.2e'
          % (np.max(np.abs(U2 - U1)) / np.max(np.abs(U1))))
    print('        （方孔对称 -> 镜像不可见；非对称场见验证 3）')
    print('    two_step("+"支) vs one_step 复场 max|diff|/max|U| = %.2e'
          % (np.max(np.abs(U2p - U1)) / np.max(np.abs(U1))))
    # 表 6.2 与网格关系
    m_t = 2.0
    d2_t = m_t * d1
    Dz1m, Dz1p = Dz / (1 - m_t), Dz / (1 + m_t)
    d1a_m = wvl * abs(Dz1m) / (N * d1)
    d1a_p = wvl * abs(Dz1p) / (N * d1)
    print('    m=2：表 6.2 核对：Dz1-/Dz = %.1f（书 -1），Dz2-/Dz = %.1f（书 2）；'
          'Dz1+/Dz = %.4f（书 1/3），Dz2+/Dz = %.4f（书 2/3）'
          % (Dz1m / Dz, (Dz - Dz1m) / Dz, Dz1p / Dz, (Dz - Dz1p) / Dz))
    print('    式 (6.19)-(6.23)：delta1a(-) = %.2f um = lambda|Dz1|/(N*d1)；'
          'delta2 复得 = lambda|Dz2|/(N*delta1a) = %.2f um = m*delta1'
          % (d1a_m * 1e6, wvl * abs(Dz - Dz1m) / (N * d1a_m) * 1e6))
    print('    式 (6.55)/(6.56)：1/(1-m) = %.6f vs d1/(d1-d2) = %.6f；'
          'm/(1-m) = %.6f vs d2/(d1-d2) = %.6f'
          % (1 / (1 - m_t), d1 / (d1 - d2_t), m_t / (1 - m_t), d2_t / (d1 - d2_t)))

    print('===== 验证 3：非对称场（m=2.5）：镜像效应与离散路径 =====')
    Na, L1a = 512, 8e-3
    d1a_ = L1a / Na
    m_a = 2.5
    d2a_ = m_a * d1a_
    Dza = 0.7
    na = np.arange(Na) - Na // 2
    xa, ya = np.meshgrid(na * d1a_, na * d1a_)
    ap_asym = (rect((xa - 0.8e-3) / 1.2e-3) * rect(ya / 2.0e-3)
               + 0.6 * circle(xa - 1.8e-3, ya + 1.2e-3, 0.5e-3))
    _, _, U_asym = ang_spec_prop(ap_asym, wvl, d1a_, d2a_, Dza)
    _, _, U_mc = ang_spec_mbranch(ap_asym, wvl, d1a_, d2a_, Dza, corrected=True)
    _, _, U_mb = ang_spec_mbranch(ap_asym, wvl, d1a_, d2a_, Dza, corrected=False)
    _, _, U_tm = two_step_prop(ap_asym, wvl, d1a_, d2a_, Dza)
    _, _, U_tp = two_step_prop_plus(ap_asym, wvl, d1a_, d2a_, Dza)
    _, _, U_tf = two_step_prop_fixed(ap_asym, wvl, d1a_, d2a_, Dza)
    pkA2 = np.max(np.abs(U_asym))
    print('    [镜像勘误 b] 书印"-"支 vs ang_spec：%.2e；'
          'vs 镜像(ang_spec) 的 |U| 差：%.2e —— 输出确为镜像'
          % (np.max(np.abs(U_tm - U_asym)) / pkA2,
             np.max(np.abs(np.abs(U_tm) - np.abs(rollflip(U_asym)))) / pkA2))
    print('    修正"-"支（翻转中间场）vs ang_spec：%.2e —— 逐位一致'
          % (np.max(np.abs(U_tf - U_asym)) / pkA2))
    print('    [离散路径] ang_spec vs (6.66)修正链：%.2e；" +"支 vs ang_spec：%.2e'
          % (np.max(np.abs(U_asym - U_mc)) / pkA2,
             np.max(np.abs(U_tp - U_asym)) / pkA2))
    print('    同一连续场、不同离散路径；随 m 变化的链间差：')
    for m_scan in (1.5, 2.0, 2.5, 3.0, 4.0):
        _, _, A_ = ang_spec_prop(ap_asym, wvl, d1a_, m_scan * d1a_, Dza)
        _, _, C_ = ang_spec_mbranch(ap_asym, wvl, d1a_, m_scan * d1a_, Dza, True)
        print('        m=%.1f：|mcorr-ang|/pk = %.3f'
              % (m_scan, np.max(np.abs(C_ - A_)) / np.max(np.abs(A_))))
    # 细网格仲裁（网格加倍，抽取偶样本回到同一观测网格）
    Nf = 2 * Na
    d1f = L1a / Nf
    naf = np.arange(Nf) - Nf // 2
    xf, yf = np.meshgrid(naf * d1f, naf * d1f)
    ap_f = (rect((xf - 0.8e-3) / 1.2e-3) * rect(yf / 2.0e-3)
            + 0.6 * circle(xf - 1.8e-3, yf + 1.2e-3, 0.5e-3))
    _, _, A_f = ang_spec_prop(ap_f, wvl, d1f, m_a * d1f, Dza)
    _, _, C_f = ang_spec_mbranch(ap_f, wvl, d1f, m_a * d1f, Dza, True)
    _, _, P_f = two_step_prop_plus(ap_f, wvl, d1f, m_a * d1f, Dza)
    jf = 2 * (na + Na // 2)                           # 细网格偶样本 = 粗网格坐标
    print('    网格加倍仲裁（N=512 -> 1024，坐标核对 max|dx| = %.1e m）：'
          'ang_spec %.3f、mcorr %.3f、two_step+ %.3f（x 峰值）—— 同量级收敛'
          % (np.max(np.abs((jf - Nf // 2) * d1f - na * d1a_)),
             np.max(np.abs(U_asym - A_f[np.ix_(jf, jf)])) / pkA2,
             np.max(np.abs(U_mc - C_f[np.ix_(jf, jf)])) / pkA2,
             np.max(np.abs(U_tp - P_f[np.ix_(jf, jf)])) / pkA2))
    # 同网格互证（m0 = one-step 原生缩放）
    m0 = (wvl * Dza / (Na * d1a_)) / d1a_
    _, _, U_one_m0 = one_step_prop(ap_asym, wvl, d1a_, Dza)
    _, _, U_tp_m0 = two_step_prop_plus(ap_asym, wvl, d1a_, m0 * d1a_, Dza)
    _, _, U_tf_m0 = two_step_prop_fixed(ap_asym, wvl, d1a_, m0 * d1a_, Dza)
    print('    同网格互证（m0 = %.2f）：one_step vs two_step" +"支 = %.2e；'
          'vs 修正"-"支 = %.2e'
          % (m0, np.max(np.abs(U_one_m0 - U_tp_m0)) / np.max(np.abs(U_tp_m0)),
             np.max(np.abs(U_one_m0 - U_tf_m0)) / np.max(np.abs(U_tf_m0))))

    print('===== 验证 4：式 (6.66) 勘误（书印链 = -conj(真值)） =====')
    print('    修正链 vs ang_spec：%.2e；书印链 vs ang_spec：%.2e'
          % (np.max(np.abs(U_mc - U_asym)) / pkA2,
             np.max(np.abs(U_mb - U_asym)) / pkA2))
    print('    书印链 vs -conj(ang_spec)：max|diff|/max|U| = %.2e'
          '（实输入场下逐点解析恒等 —— 书印链是"反传播共轭"，不是传播算子）'
          % (np.max(np.abs(U_mb + np.conj(U_asym))) / pkA2))

    print('===== 验证 5：Listing 6.7 的 ABCD 传播 =====')
    ABCD_free = np.array([[1.0, Dz], [0.0, 1.0]])
    _, _, U_af = ang_spec_propABCD(ap, wvl, d1, d2_book, ABCD_free)
    _, _, U_as = ang_spec_prop(ap, wvl, d1, d2_book, Dz)
    print('    f=inf（自由空间 ABCD）vs ang_spec：max|diff|/max|U| = %.2e'
          % (np.max(np.abs(U_af - U_as)) / np.max(np.abs(U_as))))
    f_lens = 2.0
    ABCD_lens = np.array([[1.0, Dz], [0.0, 1.0]]) @ \
        np.array([[1.0, 0.0], [-1.0 / f_lens, 1.0]])
    d2_lens = 1.5 * d1
    _, _, U_al = ang_spec_propABCD(ap, wvl, d1, d2_lens, ABCD_lens)
    # Collins(ABCD) == 先乘薄透镜相位 e^{-i k r1^2/(2f)} 再 Fresnel 传播（指数逐项相同）；
    # 取 ang_spec 链作参照（与 ABCD 链是同一离散算子，应到机器精度）
    _, _, U_lr2 = ang_spec_prop(ap * np.exp(-1j * k * (x1 ** 2 + y1 ** 2)
                                           / (2.0 * f_lens)),
                                wvl, d1, d2_lens, Dz)
    _, _, U_lr3 = two_step_prop_plus(ap * np.exp(-1j * k * (x1 ** 2 + y1 ** 2)
                                                / (2.0 * f_lens)),
                                     wvl, d1, d2_lens, Dz)
    print('    f=2 m 薄透镜 + Dz=1 m：ABCD 链 vs ang_spec(孔径 x e^(-i k r1^2/(2f)))：'
          'max|diff|/max|U| = %.2e（同一离散算子）'
          % (np.max(np.abs(U_al - U_lr2)) / np.max(np.abs(U_lr2))))
    print('        vs two_step"+"支（独立离散路径）参照：max|diff|/max|U| = %.2e'
          % (np.max(np.abs(U_al - U_lr3)) / np.max(np.abs(U_lr3))))

    print('===== 验证 6：二次相位采样判据 =====')
    f_nyq1 = 1.0 / (2 * d1)
    print('    源面啁啾（网格边缘 r=%.1f mm）：局部频率 %.0f /m vs Nyquist %.0f /m'
          '（比值 %.3f）'
          % (L / 2 * 1e3, (L / 2) / (wvl * Dz), f_nyq1,
             (L / 2) / (wvl * Dz) / f_nyq1))
    print('    缩放链的输入啁啾判据 |1-m|*r_edge/(lambda*Dz) <= 1/(2*d1)：'
          '非对称场 m=2.5（ang_spec 链，曲率 |1-m|）：比值 %.3f'
          % (abs(1 - m_a) * (L1a / 2) / (wvl * Dza) * 2 * d1a_))
    print('    (6.66)修正链的输入啁啾曲率 (1+m)：比值 %.3f —— 更陡，是链间差的主因'
          % ((1 + m_a) * (L1a / 2) / (wvl * Dza) * 2 * d1a_))
    print('    注：观测面 Q3 是 IFT 之后的逐点相位乘法，其后无 FT，不产生混叠；'
          'delta2 只需分辨场自身结构（近轴场集中在中部）。')

    if not outdir:
        plt.show()
        return

    # ---- 图 6-13（自绘）：四种传播器交叉验证 + 两处勘误 ----
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.6))
    axs[0, 0].imshow(ap, cmap='gray',
                     extent=[-L / 2 * 1e3, L / 2 * 1e3, -L / 2 * 1e3, L / 2 * 1e3])
    axs[0, 0].set_title('(a) 源面：2 mm 方孔（Listing 6.2，N=1024）')
    axs[0, 0].set_xlabel('x1 [mm]'); axs[0, 0].set_ylabel('y1 [mm]')
    sl = N // 2
    m_ = np.abs(x2[0, :]) <= 3e-3
    axs[0, 1].plot(x2[0, m_] * 1e3, Ian[sl, m_], 'k-', lw=1.6,
                   label='解析（式 6.1 闭式）')
    axs[0, 1].plot(x2[0, m_] * 1e3, I1[sl, m_], 'r.', ms=4,
                   label='one_step（式 6.5）')
    axs[0, 1].plot(x2[0, m_] * 1e3, I2[sl, m_], 'b.', ms=3,
                   label='two_step m=10.24（式 6.17）')
    axs[0, 1].set_title('(b) y2=0 辐照度切片：数值 vs 解析')
    axs[0, 1].set_xlabel('x2 [mm]'); axs[0, 1].set_ylabel('辐照度 [a.u.]')
    axs[0, 1].legend(fontsize=8); axs[0, 1].grid(alpha=0.3)
    slw = m_ & (Ian[sl, :] > 0.05 * np.max(Ian[sl, m_]))
    r1_ = U1[sl, slw] / Uan[sl, slw]
    ph1_ = np.angle(r1_) - np.mean(np.angle(r1_))
    axs[0, 2].plot(x2[0, slw] * 1e3, ph1_, 'b.', ms=3)
    axs[0, 2].set_title('(c) one_step/解析 相位残差（|U|>5%% 峰值处）\n'
                        'max %.1e rad —— 逐点一致' % np.max(np.abs(ph1_)))
    axs[0, 2].set_xlabel('x2 [mm]'); axs[0, 2].set_ylabel('相位差 [rad]')
    axs[0, 2].grid(alpha=0.3)
    ext = np.max(xa) * 1e3
    axs[1, 0].imshow(np.abs(U_asym), cmap='inferno',
                     extent=[-ext, ext, -ext, ext])
    axs[1, 0].set_title('(d) 非对称场（离轴矩形+圆盘）\nang_spec，m=2.5, Dz=0.7 m')
    axs[1, 0].set_xlabel('x2 [mm]'); axs[1, 0].set_ylabel('y2 [mm]')
    sa = Na // 2
    axs[1, 1].plot(xa[0, :] * 1e3, np.abs(U_asym)[sa, :] / pkA2, 'k-', lw=1.6,
                   label='ang_spec 真值 |U|')
    axs[1, 1].plot(xa[0, :] * 1e3, np.abs(U_tm)[sa, :] / pkA2, 'b--', lw=1.0,
                   label='Listing 6.3 书印"-"支 |U|（镜像场）')
    axs[1, 1].plot(xa[0, :] * 1e3,
                   np.abs(rollflip(U_tm))[sa, :] / pkA2, 'r.', ms=3,
                   label='书印输出镜像翻转后（回到真值）')
    axs[1, 1].set_title('(e) 勘误 b：书印两步"-"支输出镜像\n'
                        '（max|书印-真值|/pk = %.2f）'
                        % (np.max(np.abs(U_tm - U_asym)) / pkA2))
    axs[1, 1].set_xlabel('x2 [mm]'); axs[1, 1].set_ylabel('归一化 |U|')
    axs[1, 1].legend(fontsize=8); axs[1, 1].grid(alpha=0.3)
    axs[1, 2].plot(xa[0, :] * 1e3, np.real(U_asym)[sa, :] / pkA2, 'k-', lw=1.6,
                   label='真值 Re U')
    axs[1, 2].plot(xa[0, :] * 1e3, -np.real(np.conj(U_asym))[sa, :] / pkA2,
                   'b--', lw=1.0, label='书印 (6.66) 链 Re U（= -Re conj）')
    axs[1, 2].plot(xa[0, :] * 1e3, np.real(U_mc)[sa, :] / pkA2, 'r.', ms=3,
                   label='修正 (6.66) 链（与真值同离散路径）')
    axs[1, 2].set_title('(f) 勘误 a：书印 (6.66) 链输出 = -conj(真值)\n'
                        'max|书印+conj(真值)|/pk = %.1e'
                        % (np.max(np.abs(U_mb + np.conj(U_asym))) / pkA2))
    axs[1, 2].set_xlabel('x2 [mm]'); axs[1, 2].set_ylabel('归一化 Re U')
    axs[1, 2].legend(fontsize=8); axs[1, 2].grid(alpha=0.3)
    fig.suptitle('图 6-13（自绘仿真）工程 1：单步/两步/角谱/ABCD 交叉验证'
                 '（方孔 D=2 mm, λ=1 μm, Δz=1 m；非对称场 m=2.5, Δz=0.7 m）',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fp = os.path.join(outdir, 'fig-6-13.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
