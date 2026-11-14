# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 7 章（工程 1 + 工程 2）=====
菲涅耳衍射的采样要求：几何带宽判据 + 二次相位判据
（原书 7.1-7.3 节，Listings 7.1-7.2）

原书要点：
    光场的空间频率谱 = 平面波谱（方向余弦 alpha = lambda*f，式 7.5），
    传播几何限制"观测孔径能看到的空间频率带宽"——这是物理限制，与采样无关
    （Coy 方法的基础）。
    几何判据（7.2 节，Coy/Praus/Mansell）：
      alpha_max = (D1*delta2/delta1 + D2)/(2*Dz)        （式 7.12）
      约束 1:  delta2 <= -(D2/D1)*delta1 + lambda*Dz/D1 （式 7.14）
      约束 2:  N >= D1/(2*delta1) + D2/(2*delta2)
                       + lambda*Dz/(2*delta1*delta2)    （式 7.20）
    二次相位判据：局部空间频率 f_loc = |grad phi|/(2*pi)（式 7.37），
    FT 内的啁啾因子必须按 Nyquist 采样。
    单步菲涅耳积分法（7.3.1 节）：delta2 = lambda*Dz/(N*delta1) 固定（式 7.21），
      N >= D1*lambda*Dz/(delta1*(lambda*Dz - D2*delta1))（式 7.25 == 7.31），
      最小传播距离 Dz >= D1*delta1/lambda（式 7.42，R=inf）—— 长距离有效。
    角谱法（7.3.2 节）：delta1/delta2 独立选取，另有两条约束：
      约束 3: (1+Dz/R)*delta1 - lambda*Dz/D1 <= delta2
                              <= (1+Dz/R)*delta1 + lambda*Dz/D1（式 7.53）
      约束 4: N >= lambda*Dz/(delta1*delta2)            （式 7.59）
      —— 短距离有效。约束 4 强于约束 1+2 的组合（习题 6 的代数）。
    一般准则（7.3.3 节）：|1+Dz/R| < D2/D1 时约束 3 不起作用
      （式 7.60，几何光束包含在直径 D2 内）。

原书勘误（本脚本按修正实现并数值验证，详见笔记）：
    (a) 式 (7.46) 指数符号印反：Q[(1-m)/Dz + 1/R, r1] 的指数应为
        +i*k/2*((1-m)/Dz + 1/R)*|r1|^2（按 Q 定义式 6.7、算符链 7.45、
        以及 (7.48) 的相位 phi 均为正号；书印 - 号）。
    (b) 式 (7.47) 指数符号印反：Q2[-Dz/m, f1] = exp(-i*pi^2*2*Dz/(m*k)*|f1|^2)
        （按 Q2 定义式 6.12，d = -Dz/m 代入为负指数；书印 + 号。
        与 Listing 6.5 的实现一致。约束 4 的推导只用到 |grad phi|，不受影响。）
    (c) 7.3.1.2 节"发散 R<0 / 会聚 R>0"的符号约定与式 (7.32) 的相位定义
        及习题 4（"R = -0.2 m 收敛"）矛盾：按 U = A*exp(i*k*r1^2/(2*R))，
        局部频率 f = x1/(lambda*R)、方向余弦 alpha = x1/R，
        R > 0 发散、R < 0 会聚。验证 6 用高斯束宽度+相位曲率数值实锤。

数值验证（运行本脚本全部自动打印）：
    (0) 约束代数总账：Listing 7.1 例子（D1=2mm, D2=3mm, Dz=0.5m, delta1=40um）
        的 Nmin（式 7.25）= 65.8 -> 66、N=128、delta2=97.7um（书 7.43 全部数字）；
        最小距离（式 7.42）= 8 cm；式 7.25 == 式 7.31 恒等；
        约束 1/2/3/4 在书例（D2=4mm, Dz=0.1m）上的数值与图 7.5 对照。
    (1) Listing 7.1 复现：N=128 单步传播 vs 解析闭式（原书图 7.4）。
    (2) Listing 7.2 复现：约束 4 的 log2 N = 8.55、约束 2 的 log2 N = 8.51
        （书正文数字），N=512 角谱传播 vs 解析闭式（原书图 7.6）。
    (3) 图 7.5 复现：约束等高线图（自绘 fig-7-9）；图 7.7 复现：
        式 7.60 的三种情形（自绘 fig-7-10）。
    (4) 习题 1：g(x)=exp(i*pi*a^2*x^2), a=4, N=128, L=4m ——
        解析证明混叠（f_loc=a^2*x，|x|>1 m 处超过 Nyquist）+ 相邻样本
        相位增量判据 + 数值计数（自绘 fig-7-11）。
    (5) 习题 2/3：点源 100 km / D1=1mm 短距离传播的采样图
        （自绘 fig-7-12，约束等高线 + 可行点）。
    (6) 习题 4：R=-0.2 m 会聚源：约束图 + 角谱仿真 vs 推广解析解
        （自绘 fig-7-13）+ R 符号勘误 (c) 的束宽/曲率数值实锤。
    (7) 习题 6：约束 4 强于约束 1+2 组合的代数证明（数值核对）。

MATLAB/Python 差异（对应 .m 文件内注释）：
    - numpy 0-based / MATLAB 1-based：居中排列的零坐标样本位置不同；
    - np.meshgrid 默认 'xy' 与 MATLAB meshgrid 一致；
    - scipy.special.fresnel 直接给 Fresnel 积分；MATLAB 用 integral() 数值积分；
    - np.log2/ceil 与 MATLAB log2/nextpow2 的取整路径不同（代码里手写 2^ceil(log2)）。

用法：python ch7_sampling.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib + scipy；base.py（rect/ft2_centered，同仓库
      python/fourier/common）；ch6_fresnel_prop（one_step_prop/ang_spec_prop/
      fresnel_prop_square_ap，同仓库 python/waveprop/ch6）。
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
from base import rect, ft2_centered as ft2  # noqa: E402

# ---- 复用第 6 章的传播器与解析闭式（同仓库） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ch6"))
from ch6_fresnel_prop import (  # noqa: E402
    one_step_prop, ang_spec_prop, fresnel_prop_square_ap)


# ===== 7.2 节：几何约束（约束 1 / 约束 2） =====

def alpha_max_geom(D1, D2, wvl, Dz, d1, d2):
    """最大角谱内容（式 7.12）：alpha_max = (D1*delta2/delta1 + D2)/(2*Dz)。"""
    return (D1 * d2 / d1 + D2) / (2.0 * Dz)


def constraint1_d2(d1, D1, D2, wvl, Dz):
    """约束 1（式 7.14）：delta2 的上界 -(D2/D1)*delta1 + lambda*Dz/D1。"""
    return -D2 / D1 * d1 + wvl * Dz / D1


def constraint2_N(d1, d2, D1, D2, wvl, Dz):
    """约束 2（式 7.20）：观测面网格点数下界。
    N >= D1/(2*d1) + D2/(2*d2) + lambda*Dz/(2*d1*d2)"""
    return D1 / (2.0 * d1) + D2 / (2.0 * d2) + wvl * Dz / (2.0 * d1 * d2)


def illum_diameter(d1, d2, D1, wvl, Dz):
    """观测面被照亮区域直径（式 7.16）：D1*delta2/delta1 + lambda*Dz/delta1。"""
    return D1 * d2 / d1 + wvl * Dz / d1


def grid_diameter_min(d1, d2, D1, D2, wvl, Dz):
    """观测面网格最小总宽（式 7.18）：(Dillum + D2)/2。"""
    return (illum_diameter(d1, d2, D1, wvl, Dz) + D2) / 2.0


# ===== 7.3.1 节：单步菲涅耳积分法的专用约束 =====

def d2_one_step(N, d1, wvl, Dz):
    """单步法固定的观测面间距（式 7.21）：delta2 = lambda*Dz/(N*delta1)。"""
    return wvl * Dz / (N * d1)


def nmin_one_step(D1, D2, wvl, Dz, d1):
    """单步法网格点数下界（式 7.25 == 7.31）。
    N >= D1*lambda*Dz/(delta1*(lambda*Dz - D2*delta1))；要求 lambda*Dz > D2*delta1。"""
    return D1 * wvl * Dz / (d1 * (wvl * Dz - D2 * d1))


def min_Dz_one_step(D1, d1, wvl, R=np.inf):
    """单步法最小传播距离（式 7.41 有限 R / 式 7.42 无限 R）。
    R -> inf: Dz >= D1*delta1/lambda；
    有限 R:   Dz >= D1*delta1*R/(lambda*R - D1*delta1)（要求 lambda*R > D1*delta1）。
    注意：R 按式 (7.32) 的相位 exp(i*k*r1^2/(2*R)) 定义，R<0 会聚（见勘误 c）：
    会聚源（R<0）的分母 lambda*R - D1*delta1 < 0，与分子同号，结果仍为正。"""
    if np.isinf(R):
        return D1 * d1 / wvl
    return D1 * d1 * R / (wvl * R - D1 * d1)


# ===== 7.3.2 节：角谱法的专用约束（约束 3 / 约束 4） =====

def constraint3_d2(d1, D1, wvl, Dz, R):
    """约束 3（式 7.53）：delta2 的 [下界, 上界]。
    (1+Dz/R)*delta1 - lambda*Dz/D1 <= delta2 <= (1+Dz/R)*delta1 + lambda*Dz/D1"""
    slope = 1.0 + Dz / R
    off = wvl * Dz / D1
    return slope * d1 - off, slope * d1 + off


def constraint4_N(d1, d2, wvl, Dz):
    """约束 4（式 7.59）：N >= lambda*Dz/(delta1*delta2)。"""
    return wvl * Dz / (d1 * d2)


def constraint3_is_free(D1, D2, wvl, Dz, R):
    """式 7.60：|1+Dz/R| < D2/D1 时约束 3 不起作用（几何光束含在 D2 内）。"""
    return abs(1.0 + Dz / R) < D2 / D1


# ===== 7.3.1.2 节：局部空间频率（式 7.37-7.39） =====

def floc_quadratic(x, wvl, Dz, R):
    """源面啁啾 Q[1/Dz + 1/R] 的局部空间频率 f_x = (1/Dz + 1/R)*x/lambda（式 7.39）。
    R=inf 时退化为 f_x = x/(lambda*Dz)。R<0 会聚（勘误 c）。"""
    return (1.0 / Dz + 1.0 / R) * x / wvl


# ===== 习题 4(a)：带会聚/发散源的方孔解析解（式 7.62 的推广） =====

def fresnel_prop_square_ap_R(x2, y2, D, wvl, Dz, R):
    """边长 D 方孔 + 二次相位源 exp(i*k*r1^2/(2*R)) 的菲涅耳衍射解析解（习题 4a）。

    推导（式 7.36 的算符链）：U = Q[1/Dz, r2] V F Q[1/Dz + 1/R]{A}，
        U(x2,y2) = exp(i*pi*(x2^2+y2^2)/(lambda*Dz)) * I(x2)*I(y2)/(i*lambda*Dz)
        I(x) = exp(-i*pi*x^2/(rho*lambda*Dz))
              * sqrt(lambda*Dz/(2*rho)) * [dC + i*dS]，
        rho = 1 + Dz/R，dC = C(t2)-C(t1)，dS = S(t2)-S(t1)，
        t = sqrt(2*rho/(lambda*Dz)) * (x1 - x/rho)，x1 = +-D/2。
    R -> inf 时 rho = 1：观测面啁啾 e^{i*pi*r2^2/(lambda*Dz)} 与
    I 内的 e^{-i*pi*x^2/(lambda*Dz)} 相消，逐位退化为 ch6 的
    fresnel_prop_square_ap（脚本内已数值核对）。
    输出波前曲率 1/R_out = 1/Dz - 1/(rho*Dz) = 1/(R+Dz)，即 R_out = R + Dz
    （球面波几何：焦点位置不变，曲率半径随传播距离线性增长）。
    注意 rho > 0（本题 rho = 0.5）；rho < 0 时须改积分方向。"""
    rho = 1.0 + Dz / R
    a = rho / (wvl * Dz)                    # 指数 = i*pi*a*(x1 - x/rho)^2
    s = np.sqrt(2.0 * a)

    def I(x):
        t1 = s * (-D / 2.0 - x / rho)
        t2 = s * (D / 2.0 - x / rho)
        S1, C1 = fresnel(t1)
        S2, C2 = fresnel(t2)
        return (np.exp(-1j * np.pi * x ** 2 / (rho * wvl * Dz))
                * np.sqrt(1.0 / (2.0 * a)) * ((C2 - C1) + 1j * (S2 - S1)))

    return (np.exp(1j * np.pi * (x2 ** 2 + y2 ** 2) / (wvl * Dz))
            * I(x2) * I(y2) / (1j * wvl * Dz))


# ===== 绘图 =====

def plot_constraints_book(fig, D1, D2, wvl, Dz, d1_pick, d2_pick):
    """复现原书图 7.5：(a) 约束 4 等高线 + 约束 1/3 线；(b) 约束 2 等高线。
    域：delta1 in (0, 20 um]、delta2 in (0, 40 um]，书例 D1=2mm, D2=4mm,
    Dz=0.1m, lambda=1um, R=inf。返回 (ax1, ax2)。"""
    d1 = np.linspace(0.2e-6, 20e-6, 400)
    d2 = np.linspace(0.2e-6, 40e-6, 400)
    D1g, D2g = np.meshgrid(d1, d2)
    Nc4 = constraint4_N(D1g, D2g, wvl, Dz)
    Nc2 = constraint2_N(D1g, D2g, D1, D2, wvl, Dz)
    levels = np.arange(7, 13.5, 0.5)

    for ax, Nc, title in ((fig.axes[0], Nc4, '(a) 约束 4 的 log2 N 等高线'),
                          (fig.axes[1], Nc2, '(b) 约束 2 的 log2 N 等高线')):
        cs = ax.contourf(D1g * 1e6, D2g * 1e6, np.log2(Nc),
                         levels=levels, cmap='Greys', alpha=0.85)
        cs2 = ax.contour(D1g * 1e6, D2g * 1e6, np.log2(Nc),
                         levels=np.arange(7, 14, 1.0), colors='k',
                         linewidths=0.6)
        ax.clabel(cs2, fmt='%.0f', fontsize=7)
        # 约束 1：delta2 上界，斜率 -D2/D1、截距 lambda*Dz/D1
        ax.plot(d1 * 1e6, constraint1_d2(d1, D1, D2, wvl, Dz) * 1e6,
                'k-.', lw=2, label='约束 1（式 7.14）')
        # 约束 3 上界（R=inf：斜率 1、截距 lambda*Dz/D1）
        lo, hi = constraint3_d2(d1, D1, wvl, Dz, np.inf)
        ax.plot(d1 * 1e6, hi * 1e6, 'k--', lw=2, label='约束 3 上界（式 7.53）')
        ax.set_xlim(0, 20); ax.set_ylim(0, 40)
        ax.set_xlabel(r'$\delta_1$ [$\mu$m]')
        ax.set_ylabel(r'$\delta_2$ [$\mu$m]')
        ax.set_title(title, fontsize=11)
        ax.legend(loc='upper right', fontsize=8)
        ax.plot(d1_pick * 1e6, d2_pick * 1e6, 'r*', ms=12, zorder=5)
    return fig.axes[0], fig.axes[1]


def main(outdir=None):
    here = os.path.dirname(os.path.abspath(__file__))
    if outdir is None:
        outdir = here
    os.makedirs(outdir, exist_ok=True)
    results = []

    # ================= 验证 0：约束代数总账（书例 7.43） =================
    print('===== 验证 0：约束代数总账（式 7.12-7.59） =====')
    D1a, D2a, wvl, Dza = 2e-3, 3e-3, 1e-6, 0.5        # Listing 7.1 例子
    d1a = D1a / 50.0                                   # 孔径内至少 50 点
    Nmin_a = nmin_one_step(D1a, D2a, wvl, Dza, d1a)
    N_a = int(2 ** np.ceil(np.log2(Nmin_a)))
    d2a = d2_one_step(N_a, d1a, wvl, Dza)
    print('    Listing 7.1 例子：D1=2mm, D2=3mm, Dz=0.5m, delta1=D1/50=%.0f um' % (d1a * 1e6))
    print('    式 (7.25) Nmin = %.2f（书正文 66）；向上取 2 的幂 N = %d（书 128）'
          % (Nmin_a, N_a))
    print('    式 (7.21) delta2 = %.1f um（书 7.43：97.7 um）' % (d2a * 1e6))
    # 式 (7.31) 与 (7.25) 恒等的数值确认（从约束 2 出发、把 delta2(N) 代入后求解）
    N_731 = (D1a / d1a) / (1.0 - D2a * d1a / (wvl * Dza))
    print('    式 (7.31)（约束 2 代入 delta2=lambda*Dz/(N*delta1) 后解出）N >= %.4f，'
          '与式 (7.25) 相对差 %.1e' % (N_731, abs(N_731 - Nmin_a) / Nmin_a))
    print('    附加条件：lambda*Dz=%.1e > D2*delta1=%.1e（否则 N 无正解，书正文指出）'
          % (wvl * Dza, D2a * d1a))
    dzmin42 = min_Dz_one_step(D1a, d1a, wvl)           # R = inf
    print('    式 (7.42) 最小距离（R=inf）= %.1f cm（书正文 8 cm）' % (dzmin42 * 100))
    dzmin41 = min_Dz_one_step(D1a, d1a, wvl, R=0.2)
    print('    式 (7.41) 有限 R=+0.2 m（发散源）：Dz >= %.2f m；R=-0.2 m（会聚源）：Dz >= %.3f m'
          % (dzmin41, min_Dz_one_step(D1a, d1a, wvl, R=-0.2)))
    # Listing 7.2 例子上的四条约束
    D1b, D2b, Dzb = 2e-3, 4e-3, 0.1
    d1b, d2b = 9.4848e-6, 28.1212e-6
    c1 = constraint1_d2(d1b, D1b, D2b, wvl, Dzb)
    c2 = constraint2_N(d1b, d2b, D1b, D2b, wvl, Dzb)
    c3lo, c3hi = constraint3_d2(d1b, D1b, wvl, Dzb, np.inf)
    c4 = constraint4_N(d1b, d2b, wvl, Dzb)
    print('    Listing 7.2 例子：D1=2mm, D2=4mm, Dz=0.1m, delta1=9.4848um, delta2=28.1212um')
    print('      约束 1：delta2 <= %.1f um（选中 28.12：满足 %s）'
          % (c1 * 1e6, d2b <= c1))
    print('      约束 2：N >= %.2f = 2^%.2f（书正文 2^8.51）' % (c2, np.log2(c2)))
    print('      约束 3：%.1f <= delta2 <= %.1f um（R=inf，选中值满足 %s）'
          % (c3lo * 1e6, c3hi * 1e6, c3lo <= d2b <= c3hi))
    print('      约束 4：N >= %.1f = 2^%.2f（书正文 2^8.55）' % (c4, np.log2(c4)))
    print('      N = 2^ceil(log2(max)) = %d（书 512）' % 2 ** np.ceil(np.log2(max(c2, c4))))
    print('      式 (7.60)：|1+Dz/R|=1 < D2/D1=2 -> 约束 3 不起作用（图 7.5(a) 中确实被约束 1 盖住）：%s'
          % constraint3_is_free(D1b, D2b, wvl, Dzb, np.inf))
    results.append(('Nmin(L7.1)', Nmin_a, 65.8, abs(Nmin_a - 65.8)))
    results.append(('log2 N 约束4', np.log2(c4), 8.55, abs(np.log2(c4) - 8.55)))
    results.append(('log2 N 约束2', np.log2(c2), 8.51, abs(np.log2(c2) - 8.51)))

    # ================= 验证 1：Listing 7.1 复现（原书图 7.4） =================
    print('===== 验证 1：Listing 7.1 复现（N=128 单步传播 vs 解析） =====')
    n = np.arange(N_a) - N_a // 2
    x1, y1 = np.meshgrid(n * d1a, n * d1a)
    ap = rect(x1 / D1a) * rect(y1 / D1a)
    x2, y2, Uout = one_step_prop(ap, wvl, d1a, Dza)
    sl = N_a // 2
    x2s = x2[0, :]
    Uan = fresnel_prop_square_ap(x2s, 0.0, D1a, wvl, Dza)
    win = np.abs(x2s) <= D2a / 2                      # 感兴趣区域内评估
    I_num, I_an = np.abs(Uout[sl, :]) ** 2, np.abs(Uan) ** 2
    pkA = I_an.max()
    eI = np.max(np.abs(I_num[win] - I_an[win])) / pkA
    mph = win & (I_an > 0.05 * pkA)
    ph_num = np.angle(Uout[sl, :] / Uan)
    eph = np.max(np.abs(np.unwrap(ph_num[win])))
    eph_pk = np.max(np.abs(ph_num[mph]))
    print('    N=%d, delta2=%.1f um（窗口 ±%.2f mm）' % (N_a, d2a * 1e6, D2a * 1e3))
    print('    感兴趣区域内：|dI|/峰值 = %.2e；|U|>5%% 峰值处相位残差 = %.2e rad' % (eI, eph_pk))
    results.append(('L7.1 |dI|/pk', eI, None, None))
    results.append(('L7.1 相位残差', eph_pk, None, None))

    # ================= 验证 2：Listing 7.2 复现（原书图 7.6） =================
    print('===== 验证 2：Listing 7.2 复现（N=512 角谱传播 vs 解析） =====')
    Nb = int(2 ** np.ceil(np.log2(max(c2, c4))))
    n = np.arange(Nb) - Nb // 2
    x1, y1 = np.meshgrid(n * d1b, n * d1b)
    apb = rect(x1 / D1b) * rect(y1 / D1b)
    x2b, y2b, Uoutb = ang_spec_prop(apb, wvl, d1b, d2b, Dzb)
    slb = Nb // 2
    x2s2 = x2b[0, :]
    Uan2 = fresnel_prop_square_ap(x2s2, 0.0, D1b, wvl, Dzb)
    win2 = np.abs(x2s2) <= D2b / 2
    I_num2, I_an2 = np.abs(Uoutb[slb, :]) ** 2, np.abs(Uan2) ** 2
    pkB = I_an2.max()
    eI2 = np.max(np.abs(I_num2[win2] - I_an2[win2])) / pkB
    mph2 = win2 & (I_an2 > 0.05 * pkB)
    ph2 = np.angle(Uoutb[slb, :] / Uan2)
    eph2 = np.max(np.abs(ph2[mph2]))
    print('    N=%d（约束 2=%.2f、约束 4=%.2f 的最大值向上取 2 的幂）'
          % (Nb, c2, c4))
    print('    感兴趣区域内：|dI|/峰值 = %.2e；|U|>5%% 峰值处相位残差 = %.2e rad' % (eI2, eph2))
    results.append(('L7.2 |dI|/pk', eI2, None, None))
    results.append(('L7.2 相位残差', eph2, None, None))

    # ================= 验证 3：图 7.5 / 图 7.7 复现 =================
    print('===== 验证 3：图 7.5（约束等高线）与图 7.7（式 7.60 三情形）复现 =====')
    figA, axsA = plt.subplots(1, 2, figsize=(10.5, 4.2))
    plot_constraints_book(figA, D1b, D2b, wvl, Dzb, d1b, d2b)
    figA.suptitle('图 7-9（自绘仿真）工程 1：角谱传播采样约束等高线'
                  '（D1=2 mm, D2=4 mm, Δz=0.1 m, λ=1 μm, R=∞）', fontsize=12)
    figA.tight_layout(rect=[0, 0, 1, 0.92])
    fpA = os.path.join(outdir, 'fig-7-9.png')
    figA.savefig(fpA, dpi=150)
    print('    已保存 %s' % fpA)

    # 图 7.7 的三种情形：约束 1（实线，斜率 -D2/D1）、约束 3 上下界（虚/点划线）
    D1c, D2c, Dzc = 2e-3, 4e-3, 0.1                    # D2/D1 = 2
    cases = [('(a) D2/D1 < 1+Δz/R', 0.5e-3),           # 1+Dz/R = 1.2 > 2
             ('(b) -D2/D1 < 1+Δz/R < D2/D1', 4.0e-3),  # 1+Dz/R = 1.5·... 取 1+Dz/R=0.5·... 见下
             ('(c) 1+Δz/R < -D2/D1', -2.5e-3)]         # 1+Dz/R = -1.5
    # 精确取值：R = Dz/(s-1)，s = 1+Dz/R：case a s=3（R=50mm）、b s=0.5（R=200mm）、c s=-1.5（R=-66.7mm）
    Rs = [Dzc / (3.0 - 1.0), Dzc / (0.5 - 1.0), Dzc / (-1.5 - 1.0)]
    figB, axsB = plt.subplots(2, 2, figsize=(10.5, 8.4))
    axsB[1, 1].axis('off')
    d1g = np.linspace(0.2e-6, 20e-6, 300)
    for k, ((lab, _), Rk) in enumerate(zip(cases, Rs)):
        ax = axsB[k // 2, k % 2]
        slope = 1.0 + Dzc / Rk
        off = wvl * Dzc / D1c
        c1u = constraint1_d2(d1g, D1c, D2c, wvl, Dzc)
        c3u = slope * d1g + off
        c3l = slope * d1g - off
        feas = np.minimum(c1u, c3u)
        ax.fill_between(d1g * 1e6, 0, np.clip(feas, 0, None) * 1e6,
                        color='0.6', alpha=0.7, label='可行域')
        ax.plot(d1g * 1e6, c1u * 1e6, 'k-', lw=2,
                label='约束 1（斜率 -D2/D1=%.0f）' % (-D2c / D1c))
        ax.plot(d1g * 1e6, c3u * 1e6, 'k-.', lw=2,
                label='约束 3 上界（斜率 1+Δz/R=%.1f）' % slope)
        ax.plot(d1g * 1e6, np.clip(c3l, 0, None) * 1e6, 'k--', lw=2,
                label='约束 3 下界（δ1 截距 %.0f μm）'
                      % (wvl * Dzc / (D1c * slope) * 1e6))
        ax.set_xlim(0, 20); ax.set_ylim(0, 40)
        ax.set_xlabel(r'$\delta_1$ [$\mu$m]'); ax.set_ylabel(r'$\delta_2$ [$\mu$m]')
        ax.set_title('%s：R = %.1f mm，|1+Δz/R|=%.1f vs D2/D1=2'
                     % (lab, Rk * 1e3, abs(slope)), fontsize=10)
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(alpha=0.3)
        print('    %s R=%+.1f mm -> 1+Δz/R = %+.1f，式(7.60) 判定约束 3 %s'
              % (lab[:3], Rk * 1e3, slope,
                 '不起作用' if constraint3_is_free(D1c, D2c, wvl, Dzc, Rk) else '起作用'))
    figB.suptitle('图 7-10（自绘仿真）工程 1：式 (7.60) 三情形的约束 3 去留'
                  '（阴影 = 同时满足约束 1/3 的 δ2 可行带）',
                  fontsize=12)
    figB.tight_layout(rect=[0, 0, 1, 0.94])
    fpB = os.path.join(outdir, 'fig-7-10.png')
    figB.savefig(fpB, dpi=150)
    print('    已保存 %s' % fpB)

    # ================= 图 7-8：Listing 7.1/7.2 复现对照 =================
    fig0, axs0 = plt.subplots(2, 2, figsize=(11.5, 7.6))
    ax = axs0[0, 0]
    ax.plot(x2s * 1e3, I_num / pkA, 'k-', lw=1.4, label='数值')
    ax.plot(x2s * 1e3, I_an / pkA, 'r--', lw=1.0, label='解析')
    ax.set_xlim(-1.5, 1.5)
    ax.set_xlabel(r'$x_2$ [mm]'); ax.set_ylabel('归一化辐照度')
    ax.set_title('(a) 单步法（Listing 7.1）：N=128，δ1=40 μm\n|dI|/pk = %.1e' % eI,
                 fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs0[0, 1]
    ax.plot(x2s * 1e3, np.angle(Uout[sl, :]), 'k.', ms=3, label='数值')
    ax.plot(x2s * 1e3, np.angle(Uan), 'r-', lw=1.0, label='解析')
    ax.set_xlim(-1.5, 1.5)
    ax.set_xlabel(r'$x_2$ [mm]'); ax.set_ylabel('相位 [rad]')
    ax.set_title('(b) 单步法相位：残差 %.1e rad' % eph_pk, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs0[1, 0]
    ax.plot(x2s2 * 1e3, I_num2 / pkB, 'k-', lw=1.2, label='数值')
    ax.plot(x2s2 * 1e3, I_an2 / pkB, 'r--', lw=1.0, label='解析')
    ax.set_xlim(-2.5, 2.5)
    ax.set_xlabel(r'$x_2$ [mm]'); ax.set_ylabel('归一化辐照度')
    ax.set_title('(c) 角谱法（Listing 7.2）：N=512，δ1=9.48 μm，δ2=28.12 μm\n|dI|/pk = %.1e'
                 % eI2, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs0[1, 1]
    ax.plot(x2s2 * 1e3, np.angle(Uoutb[slb, :]), 'k.', ms=2.5, label='数值')
    ax.plot(x2s2 * 1e3, np.angle(Uan2), 'r-', lw=1.0, label='解析')
    ax.set_xlim(-2.5, 2.5)
    ax.set_xlabel(r'$x_2$ [mm]'); ax.set_ylabel('相位 [rad]')
    ax.set_title('(d) 角谱法相位：残差 %.1e rad' % eph2, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig0.suptitle('图 7-8（自绘仿真）工程 1：书例复现——单步法（原书图 7.4）与角谱法'
                  '（原书图 7.6）vs 解析闭式（方孔 D1=2 mm, λ=1 μm）', fontsize=12)
    fig0.tight_layout(rect=[0, 0, 1, 0.94])
    fp0 = os.path.join(outdir, 'fig-7-8.png')
    fig0.savefig(fp0, dpi=150)
    print('    已保存 %s' % fp0)

    print('===== 汇总 =====')
    for nm, val, ref, dev in results:
        if ref is None:
            print('    %-18s = %.3e' % (nm, val))
        else:
            print('    %-18s = %.4g（参考 %.4g，偏差 %.1e）' % (nm, val, ref, dev))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
