# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 8 章（工程 1）：用部分传播放宽采样约束 =====
原书第 8 章 8.1-8.4 节（吸收边界、两次/任意次部分传播、采样约束），
书 Listing 8.1 / 8.2 的 Python 等价实现 + 数值验证。

原书要点：
    角谱法只适合短距离：约束 4（N >= lambda*Dz/(delta1*delta2)，式 7.59）
    随 Dz 线性变紧，Dz 大时把 N 逼到离谱（书举 N > 4096）。约束 4 只与
    传播方法有关、与几何无关；固定 N、delta1、delta2、lambda 而缩短 Dz，
    约束 4 只会更松 —— 这就是把全长传播拆成多次"部分传播"的依据。
    解法（本章）：把 Dz 拆成 n-1 段（n 个平面），每段用角谱法，并在每个
    平面乘吸收边界（absorbing boundary）压掉网格边缘能量，抑制 wrap-around：
      - 吸收边界（8.1 节）：super-Gaussian（式 8.1）或 Tukey/余弦锥削窗
        （式 8.2），中心 ~1、边缘 ~0；Hamming/Bartlett 中心也被压低，不适用。
      - 网格间距沿 z 线性插值（式 8.8）：delta_i = (1-alpha_i)*delta1 + alpha_i*deltan，
        alpha_i = z_i/Dz（相似三角形，图 8.3）。
      - 中平面两个二次相位恰好互逆抵消（式 8.15：(m2-1)/Dz2 = (m1-1)/(m1*Dz1)），
        得到干净的连乘算符链（式 8.16-8.18）。
      - 采样约束（8.4 节）：约束 1/2 是几何的，原样；约束 3 与分段数无关
        （式 8.22 == 式 7.53，alpha_2 在推导中严格消去，习题 2）；约束 4 按段
        改写为 Delta_z_i <= min(delta1,deltan)^2 * N / lambda（式 8.24），
        由此定单段最大距离 Delta_z_max 与最少平面数 n = ceil(Dz/Dzmax) + 1。
    书例（Listing 8.2，终稿参数）：方孔平面波，D1=2mm、lambda=1um、Dz=2m、
    N=128、delta1=D1/30=66.7um、deltan=D2/30=200um、5 次部分传播（每段 0.4m）。

原书笔误/不一致（数值与图面取证，详见笔记对应 warning）：
    (a) 式 (8.2) Tukey 窗第一分支书印 "1, r >= alpha*L/2"：不等号印反
        （平顶应在中心，与图 8.1 矛盾），且第二分支 r/L 与 alpha*N/2、N/2
        单位混用（L 与 N 混写）。按图 8.1 与标准 Tukey 定义实现：
        平顶 r <= alpha*L/2，锥形 alpha*L/2 <= r <= L/2（L = N*delta 全宽）。
    (b) 书例参数三处不一致：
        - 正文写 "delta_n <= 133 um"（由 D2=4mm/30 来），Listing 8.2 写
          D2=6mm、deltan=D2/30=200um；
        - 图 8.5 整张按旧稿参数绘制：约束 1 虚线截距 500um = lambda*1m/D1
          （Dz=1m，不是正文的 2m）、斜率 -2 = -D2/D1（D2=4mm，不是 6mm），
          白叉实测 (66.0, 132.7) um —— 白叉处旧稿口径 log2N = 6.43，
          N >= 87 -> 取 2 的幂 128，与正文 "at least N = 2^7 = 128" 自洽；
        - 图 8.6 数值标记间距 ~0.13 mm（+-2mm 内约 30 个点），对应
          delta_n = 133um（正文口径），不是 Listing 8.2 的 200um。
        数值自洽的一套是 Listing 8.2（约束 2 需 N>=105、单次约束 4 需
        N>=150 均在 N=128 下成立/成立不了——分段确实必要）；正文 133um
        与 Dz=2m 组合则违反约束 2（N>=143）。以 Listing 8.2 为准复现。
    (c) "至少 n = ceil(Dz/Dzmax)+1 次部分传播"的计数口径：按表 8.2 记号 n
        是平面数（n-1 次传播），书例最少 4 次传播（5 个平面，每段 0.5m
        <= 0.569m）；Listing 8.2 实际取 5 次传播（每段 0.4m），比下限多
        一段（"可以多用部分传播"），无害。
    (d) 式 (8.16)（两次分段）只在中平面放 A[r2]，而式 (8.18)/Listing 8.1
        在每个平面（含最终观察面）都乘吸收边界——ROI 内 sg~1，无实际影响；
        代码按式 (8.18)/Listing 8.1 实现。
    (e) 表 8.2 把 alpha_i = z_i/Dz 描述成 "fractional distance from plane 1
        to plane i+1"：按表 8.3（alpha_1 = z_1/Dz = 0 对应平面 1）与
        Listing 8.1，alpha_i 对应平面 i，表格措辞差一。

数值验证（运行本脚本全部自动打印，图存 out_*.png）：
    (0) 式 (8.8) 网格插值 / 式 (8.15) 中平面相位抵消 / 式 (8.21)->(8.22)
        的 alpha 消去（习题 2 数值版，1000 组随机参数判定一致）。
    (1) 退化一致性：n=1（单段、无吸收边界）与第 6 章 ang_spec_prop 逐位一致。
    (2) Listing 8.2 复现（终稿参数）：N=128、5 段传播 vs 方孔解析闭式
        （原书图 8.6）；并按正文 133um 口径复算同题，量化约束 2 被突破的后果。
    (3) 采样约束总账 + 勘误 b 取证汇总（白叉坐标/虚线截距/斜率/标记间距）
        + 式 (8.25) 的 Delta_z_max 与最少平面数。
    (4) wrap-around 对比：Dz=8m 单次传播（违反约束 4）vs 分段 20 次 +
        吸收边界（out_fig7_wraparound.png -> 笔记 fig-8-7）。
    (5) 各平面场演化与边缘能量吸收（out_fig8_planes.png -> fig-8-8）。
    (6) 约束图复现与旧稿/终稿参数对照（out_fig9_constraints.png -> fig-8-9）。

与 MATLAB 版（matlab/waveprop/ch8/ch8_partial_prop.m）的差异：
    - numpy 0-based（零坐标样本在 N//2）/ MATLAB 1-based（N/2+1 行）；
    - meshgrid 行为一致；fft 归一化：ft2 乘 delta^2、ift2 乘 (N*deltaf)^2
      （书 Listing 2.5/2.6）；
    - Python 版 ang_spec_multi_prop_vac 增加 sg/record 可选参数（书版等价于默认）；
    - MATLAB 全 double 无整数溢出；numpy 侧 nsq**8 必须转 float
      （Windows numpy 默认 int32，8192**8 会回绕，见 super_gaussian_absorber）。

运行：C:/Anocanda/envs/machinelearn/python.exe -I ch8_partial_prop.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def plain_log(axis):
    """log 轴用纯文本刻度（避免 mathtext 的 U+2212 负号被 SimHei 顶成缺字框）。"""
    axis.set_major_formatter(FuncFormatter(lambda v, p: '%g' % v))
    axis.set_minor_formatter(FuncFormatter(lambda v, p: ''))


# ---- 复用公共基础模块与第 6 章传播器（同仓库） ----
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'fourier', 'common'))
from base import rect, ft2_centered as ft2          # noqa: E402
sys.path.insert(0, os.path.join(_HERE, '..', 'ch6'))
from ch6_fresnel_prop import ang_spec_prop, fresnel_prop_square_ap  # noqa: E402


def ift2(G, delta_f):
    """原书 Listing 2.6：定标正确的二维逆 DFT（乘 (N*delta_f)^2）。"""
    N = G.shape[0]
    return np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(G))) * (N * delta_f) ** 2


# ===== 8.1 节：吸收边界（式 8.1 / 8.2） =====

def super_gaussian_absorber(N, w_frac=0.47, order=16):
    """super-Gaussian 吸收边界（书 Listing 8.1 第 10-12 行的函数化）。

    书代码：sg = exp(-nsq.^8 / w^16)，nsq = nx^2+ny^2（索引单位）、w = 0.47*N。
    即式 (8.1) 的 n=16、sigma = w*delta = 0.47*N*delta = 0.47*L（L 为网格全宽），
    与图 8.1 的 sigma=0.45L 同量级。
    注意：nsq**8 必须用浮点算——Windows numpy 默认 int32，8192**8 会整数溢出
    回绕、sg 处处 ~1（MATLAB 是 double 无此坑）。
    """
    n = (np.arange(N) - N // 2).astype(float)
    nx, ny = np.meshgrid(n, n)
    nsq = nx ** 2 + ny ** 2
    w = w_frac * N
    return np.exp(-nsq ** (order // 2) / w ** order)


def tukey_absorber(N, alpha):
    """Tukey（余弦锥削）吸收边界，式 (8.2) 的修正实现（勘误 a）。

    书印两分支为 "1, r >= alpha*L/2" 与 "alpha*N/2 <= r/L <= N/2"：
    第一分支不等号印反（否则平顶在边缘、与图 8.1 矛盾），且 r 与 r/L、
    L 与 N 单位混用。按图 8.1 与标准 Tukey 定义实现（L = N*delta 全宽，
    r 取索引半径 * delta，这里直接用索引半径，delta 两边约掉）：
        平顶：r <= alpha*L/2           （索引半径 <= alpha*N/2）
        锥形：alpha*L/2 <= r <= L/2    （索引半径 alpha*N/2 -> N/2，余弦降到 0）
    方格四角（索引半径 > N/2）超出书公式定义域，钳为 0。
    """
    n = np.arange(N) - N // 2
    nx, ny = np.meshgrid(n, n)
    R = np.sqrt(nx ** 2 + ny ** 2)          # 索引半径（单位：格点）
    half = N / 2.0
    a_half = alpha * half
    g = np.ones((N, N))
    band = (R > a_half) & (R <= half)
    g[band] = 0.5 * (1.0 + np.cos(np.pi * (R[band] - a_half) / ((1.0 - alpha) * half)))
    g[R > half] = 0.0
    return g


# ===== 8.2 节：两分段的几何关系（式 8.3-8.8、8.10-8.15） =====

def delta_middle(d1, d3, alpha):
    """中间平面网格间距（式 8.8）：delta2 = (1-alpha)*delta1 + alpha*delta3。"""
    return (1.0 - alpha) * d1 + alpha * d3


def middle_chirp_cancels(d1, d3, Dz1, Dz2):
    """式 (8.15)：(m2-1)/Dz2 == (m1-1)/(m1*Dz1)（中平面两个 Q 因子互逆）。"""
    a = Dz1 / (Dz1 + Dz2)
    d2 = delta_middle(d1, d3, a)
    m1 = d2 / d1
    m2 = d3 / d2
    return (m2 - 1.0) / Dz2, (m1 - 1.0) / (m1 * Dz1)


# ===== 8.4 节：采样约束（约束 1/2/3 沿用第 7 章，约束 4 按段改写） =====

def constraint1_d2(d1, D1, D2, wvl, Dz):
    """约束 1（式 7.14）：delta_n 上界（几何，不随分段改变）。"""
    return -D2 / D1 * d1 + wvl * Dz / D1


def constraint2_N(d1, dn, D1, D2, wvl, Dz):
    """约束 2（式 7.20）：网格点数下界（几何，不随分段改变）。"""
    return D1 / (2 * d1) + D2 / (2 * dn) + wvl * Dz / (2 * d1 * dn)


def constraint3_d2(d1, D1, wvl, Dz, R):
    """约束 3（式 8.22 == 式 7.53）：与分段数无关（alpha 消去，习题 2）。"""
    slope = 1.0 + Dz / R
    off = wvl * Dz / D1
    return slope * d1 - off, slope * d1 + off


def dz_max_partial(d1, dn, N, wvl):
    """式 (8.24)：单段最大传播距离 Delta_z_max = min(d1,dn)^2 * N / lambda。"""
    return min(d1, dn) ** 2 * N / wvl


def n_planes_min(Dz, dzmax):
    """最少平面数：n = ceil(Dz/Dzmax) + 1（n-1 次传播；计数口径见勘误 c）。"""
    return int(np.ceil(Dz / dzmax)) + 1


# ===== Listing 8.1：任意次分段传播的角谱法（书的函数原样整理 + 可选扩展） =====

def ang_spec_multi_prop_vac(Uin, wvl, delta1, deltan, z, sg=None, record=False):
    """任意次部分传播（书 Listing 8.1 的 Python 等价）。

    输入
        Uin    : U(r1)，源平面光场 [(W/m^2)^(1/2)]
        wvl    : 光波长 [m]
        delta1 : 源平面网格间距 [m]
        deltan : 观察平面网格间距 [m]
        z      : 平面位置数组 [z2, z3, ..., zn]（不含 z1=0）[m]
    输出
        xn, yn : 观察平面坐标 [m]
        Uout   : U(rn) [(W/m^2)^(1/2)]
    可选扩展（书版等价于默认值）
        sg     : 吸收边界阵列；None 时用书的 super-Gaussian（w=0.47N, n=16）
        record : True 时额外返回每个平面（乘吸收边界后）的场列表
    与书的差异：numpy 0-based（零坐标样本在 N//2）；MATLAB 版为 1-based。
    """
    N = Uin.shape[0]                       # number of grid points
    n = np.arange(N) - N // 2             # 书：nx = -N/2 : N/2-1（1-based）
    nx, ny = np.meshgrid(n, n)
    k = 2.0 * np.pi / wvl                  # optical wavevector
    # super-Gaussian absorbing boundary（书第 9-12 行）
    if sg is None:
        sg = super_gaussian_absorber(N)
    z = np.concatenate(([0.0], np.asarray(z, dtype=float)))   # 书：z = [0 z]
    n_planes = len(z)
    Delta_z = z[1:] - z[:-1]               # propagation distances
    alpha = z / z[-1]                      # fractional distances（alpha_i = z_i/Dz）
    delta = (1.0 - alpha) * delta1 + alpha * deltan           # 式 (8.8)：各平面网格间距
    m = delta[1:] / delta[:-1]             # scaling factors m_i
    x1 = nx * delta[0]
    y1 = ny * delta[0]
    r1sq = x1 ** 2 + y1 ** 2

    Q1 = np.exp(1j * k / 2.0 * (1.0 - m[0]) / Delta_z[0] * r1sq)
    Uin = Uin * Q1
    planes = []
    for idx in range(n_planes - 1):        # 书：for idx = 1 : n-1
        # spatial frequencies (of i-th plane)
        deltaf = 1.0 / (N * delta[idx])
        fX = nx * deltaf
        fY = ny * deltaf
        fsq = fX ** 2 + fY ** 2
        Z = Delta_z[idx]                   # propagation distance
        # quadratic phase factor（式 8.18 的 Q2[-Delta_z_i/m_i]，负号已对照原书）
        Q2 = np.exp(-1j * np.pi ** 2 * 2.0 * Z / m[idx] / k * fsq)
        # compute the propagated field（含 A[r_{i+1}] 吸收边界，按式 8.18）
        Uin = sg * ift2(Q2 * ft2(Uin / m[idx], delta[idx]), deltaf)
        if record:
            planes.append(Uin.copy())
    # observation-plane coordinates
    xn = nx * delta[-1]
    yn = ny * delta[-1]
    rnsq = xn ** 2 + yn ** 2
    # 注意：循环结束后 Z 恰为最后一段 Delta_z(n-1)（书代码同款写法）
    Q3 = np.exp(1j * k / 2.0 * (m[-1] - 1.0) / (m[-1] * Z) * rnsq)
    Uout = Q3 * Uin
    if record:
        return xn, yn, Uout, planes
    return xn, yn, Uout


# ===== 主流程：验证 0-6 =====

def main(outdir=None):
    outdir = outdir or os.path.dirname(os.path.abspath(__file__))

    # ================= 验证 0：几何关系与 alpha 消去（习题 2 数值版） =================
    print('===== 验证 0：式 (8.8) 网格插值 / 式 (8.15) 相位抵消 / 式 (8.21)->(8.22) =====')
    d1, d3 = 66.7e-6, 200e-6
    a1 = 0.3
    d2 = delta_middle(d1, d3, a1)
    print('    式 (8.8)：delta1=66.7um, delta3=200um, alpha=0.3 -> delta2 = %.4f um'
          % (d2 * 1e6))
    lhs, rhs = middle_chirp_cancels(d1, d3, 0.6, 1.4)
    print('    式 (8.15)：(m2-1)/Dz2 = %.6e 与 (m1-1)/(m1*Dz1) = %.6e，相对差 %.1e'
          % (lhs, rhs, abs(lhs - rhs) / abs(rhs)))
    # 式 (8.21)->(8.22)：约束 3 的 alpha_2 消去（随机参数验证）
    rng = np.random.default_rng(8)
    mism = 0
    for _ in range(1000):
        d1r = 10 ** rng.uniform(-5, -3)
        dnr = 10 ** rng.uniform(-5, -3)
        Dzr = 10 ** rng.uniform(-1, 5)
        Rr = 10 ** rng.uniform(-1, 5) * rng.choice([-1.0, 1.0])
        wvlr = 10 ** rng.uniform(-7, -5)
        D1r = 10 ** rng.uniform(-3, -1)
        a2 = rng.uniform(0.05, 0.95)
        Dz1 = a2 * Dzr
        d2r = (1 - a2) * d1r + a2 * dnr
        lo2 = (1 + Dz1 / Rr) * d1r - wvlr * Dz1 / D1r      # 式 (8.20) 对 delta2
        hi2 = (1 + Dz1 / Rr) * d1r + wvlr * Dz1 / D1r
        ok2 = lo2 <= d2r <= hi2
        lon = (1 + Dzr / Rr) * d1r - wvlr * Dzr / D1r      # 式 (8.22) 对 deltan
        hin = (1 + Dzr / Rr) * d1r + wvlr * Dzr / D1r
        okn = lon <= dnr <= hin
        if ok2 != okn:
            mism += 1
    print('    1000 组随机参数：式 (8.20)（对 delta2）与式 (8.22)（对 deltan）'
          '判定不一致的组数 = %d —— alpha_2 严格消去，约束 3 与分段无关' % mism)

    # ================= 验证 1：n=1 退化一致性（与第 6 章 ang_spec_prop） =================
    print('===== 验证 1：n=1、无吸收边界 -> 与第 6 章 ang_spec_prop 逐位一致 =====')
    N = 128
    d1 = d3 = 100e-6
    Dz = 0.5
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    apt = rect(x1 / 2e-3) * rect(y1 / 2e-3)
    _, _, U_a = ang_spec_multi_prop_vac(apt, 1e-6, d1, d3, [Dz], sg=np.ones((N, N)))
    _, _, U_b = ang_spec_prop(apt, 1e-6, d1, d3, Dz)
    print('    最大复数差 = %.1e（应 ~1e-13 浮点位级）' % np.max(np.abs(U_a - U_b)))

    # ================= 验证 2：Listing 8.2 书例复现（原书图 8.6） =================
    print('===== 验证 2：Listing 8.2 复现（方孔平面波，N=128，5 段传播） =====')
    D1 = 2e-3
    wvl = 1e-6
    z_tot = 2.0
    delta1 = D1 / 30.0                     # 30 点跨源孔径 -> 66.67 um
    deltan = 200e-6                        # Listing 8.2：D2=6mm, deltan=D2/30（勘误 b：正文 133um）
    N = 128
    nseg = 5                               # 书：number of partial propagations
    z = (np.arange(1, nseg + 1)) * z_tot / nseg     # 书：z = (1:n)*z/n
    x1, y1 = np.meshgrid((np.arange(N) - N // 2) * delta1,
                         (np.arange(N) - N // 2) * delta1)
    ap = rect(x1 / D1) * rect(y1 / D1)
    xn, yn, Uout, planes = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z,
                                                   record=True)
    # 解析对照：y=0 切片（numpy 0-based：第 N//2 行；MATLAB 书代码为 N/2+1 行）
    sl = N // 2
    Uan = fresnel_prop_square_ap(xn[sl, :], 0.0, D1, wvl, z_tot)
    win2mm = np.abs(xn[sl, :]) <= 2e-3          # 图 8.6 绘制的 +-2mm 范围
    win_roi = np.abs(xn[sl, :]) <= 3e-3         # 终稿 ROI（D2=6mm）
    I_num, I_an = np.abs(Uout[sl, :]) ** 2, np.abs(Uan) ** 2
    pk = I_an.max()
    eI = np.max(np.abs(I_num[win_roi] - I_an[win_roi])) / pk
    mph = win_roi & (I_an > 0.05 * pk)
    eph = np.max(np.abs(np.angle(Uout[sl, :][mph] / Uan[mph])))
    print('    delta1=%.2f um, deltan=%.2f um, N=%d, 段长 %.2f m x %d'
          % (delta1 * 1e6, deltan * 1e6, N, z[0], nseg))
    mzero = (xn[sl, :] > 0.5e-3) & (xn[sl, :] < 1.5e-3)
    xzero = xn[sl, :][mzero][np.argmin(I_an[mzero])]
    print('    解析峰值辐照度 = %.2f（原书图 8.6 峰值 2.57）；第一暗环 x = %.2f mm'
          '（lambda*Dz/D1 = 1 mm 量级）' % (pk, xzero * 1e3))
    print('    ROI(+-%.0f mm)内：|dI|/峰值 = %.2e；相位残差(|U|>5%%峰值) = %.2e rad'
          % (3.0, eI, eph))
    # 正文口径（delta_n=133um）同样跑一遍：约束 2 被突破的后果
    xn1, yn1, Uout133 = ang_spec_multi_prop_vac(ap, wvl, delta1, 133e-6, z)
    Uan133 = fresnel_prop_square_ap(xn1[sl, :], 0.0, D1, wvl, z_tot)
    I_num133, I_an133 = np.abs(Uout133[sl, :]) ** 2, np.abs(Uan133) ** 2
    winb = np.abs(xn1[sl, :]) <= 2e-3
    pkb = I_an133.max()
    eI133 = np.max(np.abs(I_num133[winb] - I_an133[winb])) / pkb
    mphb = winb & (I_an133 > 0.05 * pkb)
    eph133 = np.max(np.abs(np.angle(Uout133[sl, :][mphb] / Uan133[mphb])))
    print('    [正文口径 deltan=133um] +-2mm 内：|dI|/峰值 = %.2e；相位残差 = %.2e rad'
          '（约束 2 被突破 143>128 仍大致可用——卷绕污染温和）' % (eI133, eph133))

    # ================= 验证 3：采样约束总账 + 勘误 b 取证汇总 =================
    print('===== 验证 3：书例采样约束 + 勘误 b（图 8.5 旧稿 / 正文 133um / Listing 200um）=====')
    D2 = 6e-3
    c1 = constraint1_d2(delta1, D1, D2, wvl, z_tot)
    c2_200 = constraint2_N(delta1, 200e-6, D1, D2, wvl, z_tot)
    c2_133_2m = constraint2_N(delta1, 133e-6, D1, 4e-3, wvl, z_tot)
    c2_133_1m = constraint2_N(delta1, 133e-6, D1, 4e-3, wvl, 1.0)
    lo3, hi3 = constraint3_d2(delta1, D1, wvl, z_tot, np.inf)
    print('    约束 1（D2=6mm, Dz=2m）：delta_n <= %.0f um（200um 满足）' % (c1 * 1e6))
    print('    约束 2（deltan=200um, Dz=2m）：N >= %.1f = 2^%.2f -> N=128 成立'
          % (c2_200, np.log2(c2_200)))
    print('    约束 2（deltan=133um, D2=4mm, Dz=2m，正文口径）：N >= %.1f -> N=128 不成立'
          % c2_133_2m)
    print('    约束 2（deltan=133um, D2=4mm, Dz=1m，图 8.5 旧稿口径）：N >= %.1f'
          ' = 2^%.2f -> N=128 成立（与正文"at least N=2^7"自洽）'
          % (c2_133_1m, np.log2(c2_133_1m)))
    print('    约束 3（R=inf）：%.0f <= delta_n <= %.0f um（|1+Dz/R|=1 < D2/D1=3，免检）'
          % (lo3 * 1e6, hi3 * 1e6))
    print('    约束 4（单次，deltan=200um）：N >= lambda*Dz/(d1*dn) = %.1f > 128 -> 必须分段'
          % (wvl * z_tot / (delta1 * deltan)))
    # 图 8.5 取证（PDF 180dpi 像素实测，刻度用原书刻度标签标定）：
    print('    [图 8.5 取证] 白叉实测 (66.0, 132.7) um（对应 30 点跨 D1 与 D2=4mm）；')
    print('      约束 1 虚线：delta_n 截距 ~500um = lambda*1m/D1（Dz=1m 旧稿），')
    print('        delta_n=0 交于 delta_1~250um -> 斜率 -2 = -D2/D1（D2=4mm 旧稿）；')
    print('      白叉处 log2N：旧稿口径 %.2f（N>=87->128）、正文混合口径 %.2f（N>=143->256）'
          % (np.log2(c2_133_1m), np.log2(c2_133_2m)))
    dzmax = dz_max_partial(delta1, deltan, N, wvl)
    npmin = n_planes_min(z_tot, dzmax)
    print('    式 (8.25)：Delta_z_max = min^2*N/lambda = %.4f m（书印 0.567，尾数舍入）；'
          '最少平面数 n = ceil(2/%.4f)+1 = %d（%d 次传播）'
          % (dzmax, dzmax, npmin, npmin - 1))
    print('    Listing 8.2 取 5 次传播（每段 0.4 m）—— 比下限多一段，宁多勿少')

    # ================= 验证 4：wrap-around 对比（自绘 fig-8-7） =================
    # 书例 Dz=2m、N=128 下"单次违反约束 4（需 150）"的代价对光滑方孔谱几乎
    # 不可见（sinc 谱只占 f 网格的 ~3%，边缘能量趋零）——约束违反的可见后果
    # 取决于场的谱内容。wrap-around 真正发作的场景是长距离：Dz=8m 时衍射
    # 展宽 lambda*Dz/delta1 = 120 mm 远超观察网格 25.6 mm，场大量撞上网格
    # 边缘卷绕回中心（原书动机："wrap-around creeps into the center from the
    # edge"）。对比：单次传播（无吸收边界）vs 分段 20 次（每段 0.4 m）+
    # super-Gaussian 吸收边界。
    print('===== 验证 4：长距传播的 wrap-around：单次 vs 分段 + 吸收边界 =====')
    Dz4 = 8.0
    _, _, U_single = ang_spec_prop(ap, wvl, delta1, deltan, Dz4)   # 单次、无吸收边界
    Uan4 = fresnel_prop_square_ap(xn[sl, :], 0.0, D1, wvl, Dz4)
    I_an4 = np.abs(Uan4) ** 2
    pk4 = I_an4.max()
    I_sg = np.abs(U_single[sl, :]) ** 2
    nseg4 = 20
    z4 = (np.arange(1, nseg4 + 1)) * Dz4 / nseg4
    _, _, U_part4 = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z4)
    I_pt4 = np.abs(U_part4[sl, :]) ** 2
    eI_s4 = np.max(np.abs(I_sg[win_roi] - I_an4[win_roi])) / pk4
    eI_p4 = np.max(np.abs(I_pt4[win_roi] - I_an4[win_roi])) / pk4
    print('    [Dz=8m] 约束 2 需 N >= %.0f、单次约束 4 需 N >= %.0f（均 > 128）'
          % (constraint2_N(delta1, deltan, D1, D2, wvl, Dz4),
             wvl * Dz4 / (delta1 * deltan)))
    print('      分段 20 次（0.4 m/段）：每段约束 4 需 N >= %.0f ~ %.0f（<= 128）'
          % (wvl * 0.4 / (delta1 * delta_middle(delta1, deltan, z4[0] / Dz4)),
             wvl * 0.4 / (delta_middle(delta1, deltan, z4[-2] / Dz4) * deltan)))
    print('      单次：ROI |dI|/峰值 = %.2e；分段 20 次 + 吸收边界：%.2e —— 差 %.0f 倍'
          % (eI_s4, eI_p4, eI_s4 / max(eI_p4, 1e-30)))
    ext_mm = 1e3 * deltan * N / 2
    fig7, axs7 = plt.subplots(1, 3, figsize=(13.5, 4.4))
    for ax, U, ttl in [(axs7[0], U_single, '(a) 单次传播 Δz=8m（无吸收边界）'),
                       (axs7[1], U_part4, '(b) 分段 20 次 + 吸收边界')]:
        im = ax.imshow(np.log10(np.abs(U + 1e-12)),
                       extent=[-ext_mm, ext_mm, -ext_mm, ext_mm],
                       cmap='inferno', origin='lower')
        ax.set_title(ttl + '\nlog$_{10}$|U|', fontsize=10)
        ax.set_xlabel('$x_n$ [mm]')
        fig7.colorbar(im, ax=ax, shrink=0.85)
    axs7[0].axhspan(0.4 * ext_mm, 0.95 * ext_mm, color='w', alpha=0.22)
    axs7[0].text(-ext_mm * 0.95, 0.5 * ext_mm, 'wrap-around 卷绕伪影', color='w',
                 fontsize=8)
    axs7[2].plot(xn[sl, :] * 1e3, I_an4 / pk4, 'k-', lw=1.2, label='解析')
    axs7[2].plot(xn[sl, :] * 1e3, I_pt4 / pk4, 'r.', ms=5,
                 label='分段 20 次 + 吸收')
    axs7[2].plot(xn[sl, :] * 1e3, I_sg / pk4, 'b.', ms=4, label='单次（无吸收）')
    axs7[2].set_xlim(-3, 3)
    axs7[2].set_ylim(1e-5, 2)
    axs7[2].set_yscale('log')
    plain_log(axs7[2].yaxis)
    axs7[2].set_xlabel('$x_n$ [mm]')
    axs7[2].set_ylabel('归一化辐照度')
    axs7[2].set_title('(c) ROI 内辐照度切片（对数坐标）\n单次 |dI|/pk=%.1e，分段 %.1e'
                      % (eI_s4, eI_p4), fontsize=10)
    axs7[2].legend(fontsize=8, loc='lower left')
    axs7[2].grid(alpha=0.3, which='both')
    fig7.suptitle('图 8-7（自绘仿真）长距传播的 wrap-around 与分段 + 吸收边界的抑制'
                  '（D1=2mm, D2=6mm, Δz=8m, λ=1μm, N=128）', fontsize=11)
    fig7.tight_layout(rect=[0, 0, 1, 0.92])
    fp7 = os.path.join(outdir, 'out_fig7_wraparound.png')
    fig7.savefig(fp7, dpi=150)
    print('    已保存 %s' % fp7)

    # ================= 验证 5：各平面演化与边缘能量吸收（自绘 fig-8-8） =================
    print('===== 验证 5：各平面场演化（吸收边界逐面吸收边缘能量） =====')
    nvec = np.arange(N) - N // 2
    nxs, nys = np.meshgrid(nvec, nvec)
    edge = np.sqrt(nxs ** 2 + nys ** 2) >= 0.45 * N
    edge_frac = []
    for i, Up in enumerate(planes):
        E_all = np.sum(np.abs(Up) ** 2)
        E_edge = np.sum(np.abs(Up[edge]) ** 2)
        edge_frac.append(E_edge / E_all)
        print('    平面 %d（z=%.2f m, delta=%.1f um）：总能量 %.4e，边缘占比 %.2e'
              % (i + 2, z[i], (delta1 + (z[i] / z_tot) * (deltan - delta1)) * 1e6,
                 E_all, E_edge / E_all))
    # 对照：同管线、无吸收边界（sg=1）
    _, _, _, planes0 = ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z,
                                               sg=np.ones((N, N)), record=True)
    edge_frac0 = [np.sum(np.abs(Up[edge]) ** 2) / np.sum(np.abs(Up) ** 2)
                  for Up in planes0]
    fig8, axs8 = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for i, Up in enumerate(planes):
        axs8[0].plot(nvec * delta1 * 1e3, np.abs(Up[sl, :]),
                     lw=1.1, label='平面 %d（z=%.1f m）' % (i + 2, z[i]))
    axs8[0].plot(nvec * delta1 * 1e3, np.abs(ap[sl, :]), 'k--', lw=1.0, label='源孔径')
    axs8[0].set_yscale('log')
    plain_log(axs8[0].yaxis)
    axs8[0].set_ylim(1e-4, 2)
    axs8[0].set_xlabel('索引 $\\times\\ \\delta_1$ [mm]（各平面同一索引网格）')
    axs8[0].set_ylabel('|U|（y=0 切片）')
    axs8[0].set_title('(a) 各平面 |U| 截面：边缘被 super-Gaussian 逐面剪掉', fontsize=10)
    axs8[0].legend(fontsize=7)
    axs8[0].grid(alpha=0.3, which='both')
    axs8[1].plot(range(2, 7), edge_frac0, 'bs--', ms=5,
                 label='无吸收边界（衍射能量持续外流）')
    axs8[1].plot(range(2, 7), edge_frac, 'ko-', ms=5,
                 label='有吸收边界（每段剪掉一格）')
    axs8[1].set_xlabel('平面编号')
    axs8[1].set_ylabel('边缘能量占比（r ≥ 0.45N/2）')
    axs8[1].set_title('(b) 吸收边界把边缘能量占比压低约一个量级', fontsize=10)
    axs8[1].legend(fontsize=8)
    axs8[1].grid(alpha=0.3)
    fig8.suptitle('图 8-8（自绘仿真）分段传播各平面演化：吸收边界的作用位置'
                  '（w=0.47N、n=16 super-Gaussian）', fontsize=12)
    fig8.tight_layout(rect=[0, 0, 1, 0.93])
    fp8 = os.path.join(outdir, 'out_fig8_planes.png')
    fig8.savefig(fp8, dpi=150)
    print('    已保存 %s' % fp8)

    # ================= 图 8-9：约束图复现 + 旧稿/终稿参数对照（对照原书图 8.5） ======
    d1g = np.linspace(1e-6, 260e-6, 400)
    dng = np.linspace(1e-6, 560e-6, 400)
    D1g, Dng = np.meshgrid(d1g, dng)
    Nc2_new = constraint2_N(D1g, Dng, D1, D2, wvl, z_tot)          # 终稿口径
    Nc2_old = constraint2_N(D1g, Dng, D1, 4e-3, wvl, 1.0)          # 图 8.5 旧稿口径
    fig9, ax9 = plt.subplots(figsize=(8.2, 5.6))
    cs = ax9.contourf(D1g * 1e6, Dng * 1e6, np.log2(np.maximum(Nc2_new, 1e-6)),
                      levels=np.arange(3, 14), cmap='Greys', alpha=0.55)
    cb = fig9.colorbar(cs, ax=ax9, shrink=0.9)
    cb.set_label('log2 N（约束 2，终稿 D2=6mm, Δz=2m）')
    ax9.contour(D1g * 1e6, Dng * 1e6, np.log2(np.maximum(Nc2_old, 1e-6)),
                levels=[6, 7, 8], colors='C0', linewidths=1.0,
                linestyles='dotted')
    ax9.plot(d1g * 1e6, constraint1_d2(d1g, D1, D2, wvl, z_tot) * 1e6,
             'k--', lw=2, label='约束 1（终稿：截距 1000μm）')
    ax9.plot(d1g * 1e6, constraint1_d2(d1g, D1, 4e-3, wvl, 1.0) * 1e6,
             'C0--', lw=1.6, label='约束 1（图 8.5 实绘：D2=4mm, Δz=1m，截距 500μm）')
    ax9.plot(delta1 * 1e6, 200.0, 'r*', ms=16, mec='k', zorder=6,
             label='Listing 8.2：δn=200μm（D2=6mm/30）')
    ax9.plot(66.7, 133.0, 'w*', ms=14, mec='k', zorder=6,
             label='图 8.5 白叉实测：(66.0, 132.7)μm（D2=4mm/30，旧稿）')
    ax9.annotate('白叉处 log2N=6.43（旧稿口径）\n→ N≥87 → 取 2 的幂 128',
                 (66.7, 133.0), textcoords='offset points', xytext=(20, 30),
                 fontsize=8, arrowprops=dict(arrowstyle='->', lw=0.8))
    ax9.annotate('200μm：N≥%.0f（满足）' % c2_200, (delta1 * 1e6, 200.0),
                 textcoords='offset points', xytext=(24, -30), fontsize=8,
                 arrowprops=dict(arrowstyle='->', lw=0.8))
    ax9.set_xlim(0, 260)
    ax9.set_ylim(0, 560)
    ax9.set_xlabel('$\\delta_1$ [μm]')
    ax9.set_ylabel('$\\delta_n$ [μm]')
    ax9.set_title('图 8-9（自绘）约束图复现与勘误 b：原书图 8.5 按 D2=4mm、Δz=1m 旧稿绘制，'
                  '正文/Listing 8.2 是 D2=6mm、Δz=2m', fontsize=10)
    ax9.legend(fontsize=7.5, loc='upper right')
    fig9.tight_layout()
    fp9 = os.path.join(outdir, 'out_fig9_constraints.png')
    fig9.savefig(fp9, dpi=150)
    print('    已保存 %s' % fp9)

    results = [
        ('L8.2 分段 |dI|/pk', eI, None, None),
        ('L8.2 分段 相位残差 [rad]', eph, None, None),
        ('133um 口径 |dI|/pk', eI133, None, None),
        ('Dz=8m 单次 |dI|/pk', eI_s4, None, None),
        ('Dz=8m 分段+吸收 |dI|/pk', eI_p4, None, None),
    ]
    return results


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
