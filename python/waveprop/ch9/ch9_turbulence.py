# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 9 章（工程 1）=====
通过大气湍流的传播：分步束传播法 + 相位屏 + 采样 + 验证
（原书 9.1-9.6 节，Listings 9.1-9.8）

原书要点：
    分步束传播法（split-step BPM，式 9.1-9.3）：把"介质里的传播"拆成
      真空衍射算符 R 与折射（相位）算符 T = exp(+i*phi)，phi(r) = k∫δn dz
      —— 大气折射率接近 1，只需在真空传播器（第 8 章式 8.18）的每个
      部分传播平面上乘一块相位屏。
    Kolmogorov 湍流统计（9.2.1）：惯性子区内速度/温度/折射率结构函数
      按 r^(2/3) 增长（式 9.5/9.7/9.13）；折射率 PSD（式 9.15-9.18）；
      泰勒冻结湍流假设把空间统计变成时间统计（式 9.19）。
    Rytov 微扰（9.2.2-9.2.3）：U = U0*exp(psi)，psi = chi + i*phi（式 9.25-9.27）；
      平面波相干因子、波结构函数、r0、theta0、log 振幅方差（式 9.37-9.64）。
      平面波 Kolmogorov：D_phi(r) = 6.88 (r/r0)^(5/3)（式 9.44），
      r0,pw = [0.423 k^2 Cn2 Dz]^(-3/5)（式 9.42）。
    分层模型（9.2.4-9.2.5）：n 块薄相位屏匹配连续湍流的低阶矩（式 9.65），
      各屏有效相干直径 r0i（式 9.70），由矩阵方程（式 9.75）+ 约束优化求解。
    蒙特卡洛相位屏（9.3）：FT 法（McGlamery，式 9.76-9.79）低频欠采样 ->
      Lane 次谐波法补救（式 9.81，Listing 9.3）。
    采样约束（9.4）：Coy 的"光栅模糊"模型把湍流展宽折进约束 1/2：
      D1' = D1 + c*lambda*Dz/r0,rev，D2' = D2 + c*lambda*Dz/r0（式 9.84/9.85，
      c=2 约 97% 光、c=4 约 99%），约束 3 不变；再由
      Dz_max = min(delta1,deltan)^2*N/lambda（式 9.89）定部分传播步数。

原书勘误（本脚本按修正实现并数值验证，详见笔记）：
    (a) 式 (9.2) 印成 T = exp(-i*phi)：按薄相位屏物理（时间约定
        exp(-i2*pi*nu*t)、抛物线方程 dU/dz = i/(2k)∇²U + i*k*δn*U）与
        Listing 9.8 的实现（t = exp(+i*phz)），应为 exp(+i*phi)。
        对零均值随机屏无统计影响，对确定性 δn（热透镜等）符号会反。
    (b) 式 (9.34) 的核函数：书印 sin(kr)/(kr)（3D 风格）。把式 (9.44)
        代入算出 0.568 而非 (9.49) 的 0.49（习题 3）；对 2D 各向同性 PSD，
        正确的逆变换是 J1 核 Phi(kappa) = 1/(4*pi*kappa)∫J1(kappa*r)D'(r)*r dr，
        用它恰好得 0.4896 ≈ 0.49。
    (c) 式 (9.48) 系数：书印 2*pi^2*k^2*Dz，应为 2*pi*k^2*Dz。
        反证：与 (9.42)+(9.16) 联立应给 (9.49) 的 0.49：
        2*pi*0.033/0.423 = 0.490（对）；2*pi^2 给 1.545（错，差 pi 倍）。
        且 (9.45) 的系数 6.16 = 8*pi^2*0.033/0.423（其中 8*pi^2 = 4*pi*2*pi），
        数值验证：vK 闭式 (9.45) vs 直接 Bessel 积分（用 2*pi 因子）比值恒 1.0000。
    (d) 正文 9.5.4 说"重复 40 次"，Listing 9.8 却写 nreals = 20 —— 两处不一致；
        本脚本取 40（与正文及验证图的统计量一致）。
    (e) Listing 9.8 第 31 行续行符印成 4 个点 "...."，应为 "..."。
    (f) Listing 9.6 第 21 行 d1*d2*N/wvl 与第 22 行 zmax 重复
        （d1=d2 时 min 不起作用，两式相同），非笔误、仅冗余。

实现说明（与书代码的对应关系）：
    - str_fcn2_ft / corr2_ft：书第 3 章 Listing 3.7/3.10（ch3 PDF 缺失）。
      按 9.5.5 节的用法重写：str_fcn2_ft 按定义算无偏"逐对差分"结构函数
      D(rho) = <(phi(r+rho)-phi(r))^2>（垫零到 2N 的 FFT + 重叠样本数归一）。
      注意不能用 2(var-C(rho)) 的"方差-协方差"形式：次谐波屏含波长超过
      屏幕的大幅低频模，全局方差不等于逐对平均的均方，会在小间隔处
      高估 D 数倍（本脚本开发时实测并弃用）。corr2_ft 用相关定理
      IFT[FT(u1*mask)·conj(FT(u2*mask))] 算互相干函数，第二个返回值是
      各像素对应的空间间隔模长（供径向分箱）。
    - ft_phase_screen 的合成只取 IFT 的实部：cn = (randn+i*randn)*sqrt(PSD)*df
      的 |cn|^2 期望是式 (9.79) 的 2 倍，但 Re[.] 恰好带来 1/2 因子，
      两者相消 —— Re/Im 各自都是方差正确的独立相位屏（书 Listing 9.2 第 20 行）。
    - 随机数用 numpy 的 default_rng(9) 播种，保证结果可复现（书未播种）。

数值验证（运行本脚本全部自动打印，Python 3.11 + numpy 实跑）：
    (0) 大气参数总账：r0,pw = 7.03 cm、r0,sw = 12.7 cm（书正文）、
        sigma_chi^2,sw = 0.436（书正文）；式 (9.75) 矩阵条目抽查
        alpha^(5/3) = {0, 0.0992, 0.315, 0.619, 1}（书式 9.75 一致）。
    (1) 约束优化解屏幕 r0i：目标 ||A*X-b||^2 收敛，回代 r0,sw / Rytov 数
        与目标值相对差 < 1e-6。
    (2) 采样分析：约束 2 给 log2 N = 8.47 -> N = 512（书正文"more than 2^8.47"）；
        Dz_max = 51.2 km -> nmin = 2 平面（书正文），实际用 11 平面。
    (3) 真空仿真（Listing 9.7）：ROI 内辐照度起伏 std/mean、准直后相位 RMS
        （书图 9.7：均匀辐照度 + 平坦相位）。
    (4) 湍流仿真 40 次实现（Listing 9.8）：单实现辐照度/相位（书图 9.8）。
    (5) 验证 A：末平面 40 块屏的平均结构函数 vs 6.88(r/r0)^(5/3)（书图 9.9）。
    (6) 验证 B：观测面相干因子 vs exp[-3.44(rn/r0)^(5/3)]（书图 9.10）。
    (7) FT 屏 vs 次谐波屏（Listing 9.4 参数：2 m 屏、r0=10 cm、l0=1 cm、
        L0=100 m）：低频缺失导致 FT 屏大间隔结构函数明显低于理论，
        次谐波屏贴合（书图 9.3 的复现）。

MATLAB/Python 差异（对应 .m 文件内注释）：
    - MATLAB 1-based / numpy 0-based：零坐标样本在 N/2+1 vs N//2；
    - meshgrid 行为一致（'xy'）；
    - MATLAB randn 无播种，numpy 用 default_rng(9)；
    - fmincon(约束优化) 对应 scipy.optimize.minimize(method='L-BFGS-B')；
    - 书 Listing 9.1 的 nsq.^8/w^16 直接幂可能溢出，本实现写成语义等价的
      exp(-(nsq/w^2)^8)。
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ---- 复用《计算傅里叶光学》笔记的公共基础模块（同仓库） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import ft2_centered as ft2  # noqa: E402


# ===== 书 Listing 3.2：ift2（定标正确的二维逆 DFT） =====
def ift2(G, delta_f):
    """连续逆傅里叶变换的 DFT 近似（居中输出）：
        g(r) ~= (N*delta_f)^2 * IFFT{G}
    numpy.fft.ifft2 自带 1/N^2，乘 (N*delta_f)^2 回到连续量纲。"""
    N = G.shape[0]
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(G))) * (N * delta_f) ** 2


# ===== 书 Listing 9.1：ang_spec_multi_prop =====
def ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t=None):
    """分步束传播（书 Listing 9.1）：多平面角谱传播 + 每平面相位屏 t。

    参数
    ----
    Uin      : (N,N) 复数源场
    wvl      : 波长 [m]
    delta1   : 源面网格间距 [m]
    deltan   : 观测面网格间距 [m]（中间平面按 z 线性插值）
    z        : 中间平面位置（不含 0），长度 n-1 -> 共 n 个平面
    t        : (N,N,n) 相位屏透射率 exp(i*phi)；None 表示真空传播

    返回 (xn, yn, Uout)：观测面坐标与场（全局相位 e^{ikDz} 按书约定省略）。

    式 (9.3) 的算符链：Q1 -> 循环{IFT[Q2*FT[U/m]] * T_{i+1} * sg} -> Q3。
    MATLAB/Python：零坐标样本 MATLAB 在 N/2+1，numpy 在 N//2。
    """
    N = Uin.shape[0]                    # number of grid points
    n1 = np.arange(N) - N // 2          # MATLAB: (-N/2 : N/2-1)
    nx, ny = np.meshgrid(n1, n1)
    k = 2 * np.pi / wvl
    # super-Gaussian 吸收边界（书 Listing 9.1 第 9-12 行）
    # 书：sg = exp(-nsq.^8/w^16), nsq = nx^2+ny^2, w = 0.47*N
    # 语义等价改写（避免大数幂溢出）：sg = exp(-(nsq/w^2)^8)
    nsq = nx ** 2 + ny ** 2
    w = 0.47 * N
    sg = np.exp(-(nsq / w ** 2) ** 8)

    z = np.concatenate([[0.0], np.asarray(z, dtype=float)])   # z = [0 z]
    n = len(z)
    Delta_z = z[1:] - z[:-1]            # 各部分传播距离
    alpha = z / z[-1]                   # 归一化位置（第 8 章式 8.10）
    delta = (1 - alpha) * delta1 + alpha * deltan   # 各平面网格间距（线性插值）
    m = delta[1:] / delta[:-1]          # 相邻平面缩放比
    x1 = nx * delta[0]
    y1 = ny * delta[0]
    r1sq = x1 ** 2 + y1 ** 2
    Q1 = np.exp(1j * k / 2 * (1 - m[0]) / Delta_z[0] * r1sq)
    Uin = Uin * Q1
    if t is not None:
        Uin = Uin * t[:, :, 0]          # 源平面上的相位屏 T[z1,z2] 的一部分
    for idx in range(n - 1):
        deltaf = 1 / (N * delta[idx])   # 第 i 平面的频率采样
        fX = nx * deltaf
        fY = ny * deltaf
        fsq = fX ** 2 + fY ** 2
        Z = Delta_z[idx]                # 本步传播距离
        Q2 = np.exp(-1j * np.pi ** 2 * 2 * Z / m[idx] / k * fsq)   # 频域啁啾
        Uin = ft2(Uin / m[idx], delta[idx])          # FT + 缩放
        Uin = Q2 * Uin
        Uin = ift2(Uin, deltaf)                       # IFT
        Uin = Uin * sg                                # 吸收边界
        if t is not None:
            Uin = Uin * t[:, :, idx + 1]              # 到达平面上的相位屏
    xn = nx * delta[-1]
    yn = ny * delta[-1]
    rnsq = xn ** 2 + yn ** 2
    Q3 = np.exp(1j * k / 2 * (m[-1] - 1) / (m[-1] * Z) * rnsq)  # 输出面啁啾
    Uout = Q3 * Uin
    return xn, yn, Uout


# ===== 书 Listing 9.2：ft_phase_screen =====
def ft_phase_screen(r0, N, delta, L0, l0, rng):
    """FT 法相位屏（书 Listing 9.2，式 9.76-9.79）。

    修正 von Karman 相位 PSD（式 9.51 换算到 cycles/m，式 9.52）：
        Phi_phi(f) = 0.023 r0^(-5/3) exp(-(f/fm)^2) / (f^2+f0^2)^(11/6)
    fm = 5.92/(l0*2*pi)（内尺度频率），f0 = 1/L0（外尺度频率）。
    l0=0 / L0=inf 时退化为纯 Kolmogorov（式 9.52）。
    返回一块实值相位屏 [rad]。
    """
    del_f = 1 / (N * delta)
    f1 = (np.arange(N) - N // 2) * del_f
    fx, fy = np.meshgrid(f1, f1)
    f = np.hypot(fx, fy)
    fm = 5.92 / l0 / (2 * np.pi) if l0 > 0 else np.inf    # inner scale freq
    f0 = 1.0 / L0 if np.isfinite(L0) else 0.0             # outer scale freq
    PSD_phi = 0.023 * r0 ** (-5 / 3) * np.exp(-(f / fm) ** 2) \
        / (f ** 2 + f0 ** 2) ** (11 / 6)
    PSD_phi[N // 2, N // 2] = 0          # 零频分量置零（去直流）
    # 随机傅里叶系数：cn = (randn + i*randn) * sqrt(PSD) * del_f
    # <|cn|^2> = 2*PSD*del_f^2（是式 9.79 的 2 倍），但下面只取 IFT 的实部，
    # Re 带来的 1/2 因子恰好抵消 —— Re/Im 是两块独立的合格相位屏。
    cn = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) \
        * np.sqrt(PSD_phi) * del_f
    phz = np.real(ift2(cn, 1))           # 合成（IFT 因子 N^2 即级数求和）
    return phz


# ===== 书 Listing 9.3：ft_sh_phase_screen =====
def ft_sh_phase_screen(r0, N, delta, L0, l0, rng):
    """次谐波增广相位屏（书 Listing 9.3，Lane et al.，式 9.81）。

    phz_hi：FT 法高频屏；phz_lo：3 个次谐波网格（间距 1/(3^p L)，各 3x3
    频点）的低频屏之和。返回 (phz_lo, phz_hi)，用 phz_lo + phz_hi 得整屏。
    """
    D = N * delta
    phz_hi = ft_phase_screen(r0, N, delta, L0, l0, rng)
    n1 = np.arange(N) - N // 2
    x, y = np.meshgrid(n1 * delta, n1 * delta)
    phz_lo = np.zeros((N, N))
    for p in (1, 2, 3):
        del_f = 1 / (3 ** p * D)         # 次谐波网格间距
        f1 = np.array([-1.0, 0.0, 1.0]) * del_f
        fx, fy = np.meshgrid(f1, f1)
        f = np.hypot(fx, fy)
        fm = 5.92 / l0 / (2 * np.pi) if l0 > 0 else np.inf
        f0 = 1.0 / L0 if np.isfinite(L0) else 0.0
        PSD_phi = 0.023 * r0 ** (-5 / 3) * np.exp(-(f / fm) ** 2) \
            / (f ** 2 + f0 ** 2) ** (11 / 6)
        PSD_phi[1, 1] = 0                # 3x3 网格的中心（零频）置零
        cn = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) \
            * np.sqrt(PSD_phi) * del_f
        # SH = sum cn(ii) exp(i*2*pi*(fx(ii)*x + fy(ii)*y))，9 个频点；
        # 指数可分解为内外积（数学同书 Listing 9.3 的逐点循环，只是更快）
        SH = np.zeros((N, N), dtype=complex)
        for ii in range(9):
            ex = np.exp(2j * np.pi * fx.flat[ii] * x[0, :])
            ey = np.exp(2j * np.pi * fy.flat[ii] * y[:, 0])
            SH += cn.flat[ii] * np.outer(ey, ex)
        phz_lo = phz_lo + SH             # 累加各次谐波网格
    phz_lo = np.real(phz_lo) - np.mean(np.real(phz_lo))
    return phz_lo, phz_hi


# ===== 书 Listing 3.7 的等价实现：str_fcn2_ft =====
def str_fcn2_ft(phz, delta):
    """FT 法二维结构函数（书 ch3 Listing 3.7 的等价重写，ch3 PDF 缺失）。

    直接按定义计算无偏逐对差分（不经 var-C 相减，避免大幅低频模式下
    "全局方差 - 样本协方差" 的对消误差）：
        D(rho) = <(phi(r+rho)-phi(r))^2>
               = [S2(rho) + S2(-rho) - 2*C(rho)] / N_overlap(rho)
    其中 S2/C 用垫零到 2N 的 FFT（线性、非循环）计算。
    返回 (D 二维数组, mean_r = 各像素的空间间隔模长 [m])。
    """
    N = phz.shape[0]
    q = np.zeros((2 * N, 2 * N))
    q[:N, :N] = phz                                   # 垫零 -> 循环=线性相关
    ones = np.zeros((2 * N, 2 * N))
    ones[:N, :N] = 1.0
    Fp = np.fft.fft2(q)
    Fo = np.fft.fft2(ones)
    F2 = np.fft.fft2(q ** 2)
    # 约定：ifft2(FFT(f)*conj(FFT(g)))[n] = Σ_m f[m] g[m-n]（模 2N）
    C = np.fft.ifft2(Fp * np.conj(Fp)).real    # C[n]  = Σ_{m>=n} q[m]q[m-n]
    S2b = np.fft.ifft2(F2 * np.conj(Fo)).real  # S2b[n] = Σ_{m>=n} q[m]^2（对的第 2 元）
    S2a = np.fft.ifft2(Fo * np.conj(F2)).real  # S2a[n] = Σ_{m<N-n} q[m]^2（对的第 1 元）
    ov = np.outer(N - np.arange(N), N - np.arange(N)).astype(float)  # 重叠样本数
    ov[ov == 0] = 1
    D = (S2a[:N, :N] + S2b[:N, :N] - 2 * C[:N, :N]) / ov
    X, Y = np.meshgrid(np.arange(N) * delta, np.arange(N) * delta)
    mean_r = np.hypot(X, Y)
    return D, mean_r


# ===== 书 Listing 3.10 的等价实现：corr2_ft =====
def corr2_ft(u1, u2, mask, delta):
    """FT 法互相关（书 ch3 Listing 3.10 的等价重写，ch3 PDF 缺失）。

    相关定理：Gamma(rho) = IFT[FT(u1*mask) * conj(FT(u2*mask))]。
    返回 (mcf 二维互相干函数, mean_r = 各像素的空间间隔模长 [m])。
    """
    N = u1.shape[0]
    U1 = ft2(u1 * mask, delta)
    U2 = ft2(u2 * mask, delta)
    mcf = ift2(U1 * np.conj(U2), 1 / (N * delta))
    n1 = np.arange(N) - N // 2
    X, Y = np.meshgrid(n1 * delta, n1 * delta)
    mean_r = np.hypot(X, Y)
    return mcf, mean_r


def circle_mask(X, Y, R):
    """单位高度圆盘掩模（半径 R，书 Table 1.2 的 circ）。"""
    return np.where(np.hypot(X, Y) < R, 1.0, 0.0)


# ===== 书 Listing 9.5：几何 / 大气参数 / 屏幕优化 =====
def example_setup(verbose=True):
    """书 Listing 9.5：点源 50 km、Cn2=1e-16、11 屏的几何与大气参数。

    返回 dict（几何、r0pw/r0sw/Rytov 数、屏幕 r0i、约束优化诊断）。
    """
    D2 = 0.5          # 观测望远镜直径 [m]
    wvl = 1e-6        # 波长 [m]
    k = 2 * np.pi / wvl
    Dz = 50e3         # 传播距离 [m]
    # 用 sinc 主瓣模拟点源：DROI 直径内均匀照明（书 10-11 行）
    DROI = 4 * D2     # 观测面感兴趣区直径 [m]
    D1 = wvl * Dz / DROI   # sinc 主瓣半宽的 2 倍（零点间距）[m]
    R = Dz            # 发散球面波波前曲率 [m]

    Cn2 = 1e-16       # 常数结构参数 [m^-2/3]
    # 式 (9.42)/(9.43)：常数 Cn2 时 r0,pw 与 r0,sw
    r0pw = (0.423 * k ** 2 * Cn2 * Dz) ** (-3 / 5)
    r0sw = (0.423 * k ** 2 * Cn2 * (3 / 8) * Dz) ** (-3 / 5)
    # 式 (9.64)：sigma_chi^2,sw（书用数值积分，这里同样做）
    p = np.linspace(0, Dz, 1000)
    rytov = 0.563 * k ** (7 / 6) * np.sum(Cn2 * (1 - p / Dz) ** (5 / 6)
                                          * p ** (5 / 6) * (p[1] - p[0]))

    # ---- 相位屏 r0i 的约束优化（式 9.65-9.75，书 25-41 行）----
    nscr = 11
    A = np.zeros((2, nscr))
    alpha = np.arange(nscr) / (nscr - 1)
    A[0, :] = alpha ** (5 / 3)                       # r0,sw 行（式 9.72）
    A[1, :] = (1 - alpha) ** (5 / 6) * alpha ** (5 / 6)   # Rytov 行（式 9.74）
    b = np.array([r0sw ** (-5 / 3), rytov / 1.33 * (k / Dz) ** (5 / 6)])
    x0 = (nscr / 3 * r0sw * np.ones(nscr)) ** (-5 / 3)      # 初值
    rmax = 0.1                                       # 每屏 Rytov 贡献上限
    x2 = rmax / 1.33 * (k / Dz) ** (5 / 6) / A[1, :]        # 上界（式 9.74 反解）
    x2[A[1, :] == 0] = 50 ** (-5 / 3)                # 端点屏幕给名义上界 r0=50 m
    fun = lambda X: np.sum((A @ X - b) ** 2)
    res = minimize(fun, x0, method='L-BFGS-B',
                   bounds=[(0, xx) for xx in x2])
    X = res.x
    r0scrn = X ** (-3 / 5)
    r0scrn[~np.isfinite(r0scrn)] = 1e6               # X=0 -> r0=inf（书 44 行）
    bp = A @ X
    r0sw_chk = bp[0] ** (-3 / 5)
    rytov_chk = bp[1] * 1.33 * (Dz / k) ** (5 / 6)

    if verbose:
        print('===== 验证 0：大气参数总账（书 9.5.1 节数字）=====')
        print('    D1 = %.4f m（sinc 主瓣），D2 = %.1f m，Dz = %.0f km'
              % (D1, D2, Dz / 1e3))
        print('    r0,pw = %.3f cm   r0,sw = %.3f cm（书：12.7 cm）'
              % (r0pw * 100, r0sw * 100))
        print('    r0,sw / r0,pw = %.4f（理论 (8/3)^(3/5) = %.4f）'
              % (r0sw / r0pw, (8 / 3) ** (3 / 5)))
        print('    sigma_chi^2,sw = %.3f（书：0.436）' % rytov)
        print('    式 (9.75) 矩阵第 1 行（alpha^(5/3)，5 屏情形书值）：',
              np.round(alpha[:5] ** (5 / 3), 4), '（书：0 0.0992 0.315 0.619 1）')
        print('===== 验证 1：屏幕 r0i 约束优化（Listing 9.5）=====')
        print('    exit: %s, ||A*X-b||^2 = %.3e' % (res.success, res.fun))
        print('    r0scrn [m] =', np.array2string(r0scrn, precision=3))
        print('    回代：r0,sw = %.4f cm（目标 %.4f），Rytov = %.4f（目标 %.4f）'
              % (r0sw_chk * 100, r0sw * 100, rytov_chk, rytov))
        print('    相对差：%.2e / %.2e'
              % (abs(r0sw_chk - r0sw) / r0sw, abs(rytov_chk - rytov) / rytov))
    return dict(D2=D2, wvl=wvl, k=k, Dz=Dz, DROI=DROI, D1=D1, R=R, Cn2=Cn2,
                r0pw=r0pw, r0sw=r0sw, rytov=rytov, nscr=nscr, alpha=alpha,
                r0scrn=r0scrn, res=res, r0sw_chk=r0sw_chk, rytov_chk=rytov_chk)


# ===== 书 Listing 9.6：采样分析 =====
def sampling_analysis(S, c=2.0, d1=10e-3, dn=10e-3, N=512, verbose=True):
    """书 Listing 9.6：把湍流展宽折进约束 1/2（式 9.84-9.88）并核算约束 4。

    返回 dict（D1p/D2p、约束 1/3 的边界、约束 2 的 log2N 网格、约束 4 核算）。
    """
    wvl, Dz, D1, D2, R = S['wvl'], S['Dz'], S['D1'], S['D2'], S['R']
    r0sw = S['r0sw']
    D1p = D1 + c * wvl * Dz / r0sw          # 式 (9.84)：r0,rev = r0,sw（往反路）
    D2p = D2 + c * wvl * Dz / r0sw          # 式 (9.85)
    # 约束 1（式 9.86）：deltan 上界
    dn_max = lambda d1v: (wvl * Dz - D2p * d1v) / D1p
    # 约束 3（式 9.88）：与真空相同
    lo3 = lambda d1v: (1 + Dz / R) * d1v - wvl * Dz / D1p
    hi3 = lambda d1v: (1 + Dz / R) * d1v + wvl * Dz / D1p
    # 约束 2（式 9.87）：N >= (lambda*Dz + D1p*dn + D2p*d1)/(2*d1*dn)
    N2 = lambda d1v, dnv: (wvl * Dz + D1p * dnv + D2p * d1v) / (2 * d1v * dnv)
    # 约束 4（式 9.89/9.90）
    zmax = min(d1, dn) ** 2 * N / wvl
    nmin = int(np.ceil(Dz / zmax)) + 1
    if verbose:
        print('===== 验证 2：采样分析（Listing 9.6，c = %g）=====' % c)
        print('    D1p = %.4f m, D2p = %.4f m（c*lambda*Dz/r0 = %.3f m）'
              % (D1p, D2p, c * wvl * Dz / r0sw))
        n_need = N2(d1, dn)
        print('    选点 delta1 = delta_n = %.0f mm：约束 2 需 N >= %.1f'
              '（log2 = %.2f，书：more than 2^8.47）' % (d1 * 1e3, n_need, np.log2(n_need)))
        print('    约束 1 上界 dn_max(delta1) = %.1f mm；约束 3 区间 [%.1f, %.1f] mm'
              % (dn_max(d1) * 1e3, lo3(d1) * 1e3, hi3(d1) * 1e3))
        print('    约束 4：Dz_max = %.1f km -> nmin = %d 平面（书：最小 2 平面，'
              '实际用 11 平面/10 次传播）' % (zmax / 1e3, nmin))
    return dict(D1p=D1p, D2p=D2p, dn_max=dn_max, lo3=lo3, hi3=hi3, N2=N2,
                zmax=zmax, nmin=nmin, c=c, d1=d1, dn=dn, N=N)


# ===== 书 Listing 9.7：真空仿真 =====
def vacuum_simulation(S, d1, dn, N, nplanes):
    """书 Listing 9.7：模型点源的真空传播 + 准直。

    返回 (xn, yn, Uvac, pt, sg, z)。
    """
    wvl, Dz, D1, R = S['wvl'], S['Dz'], S['D1'], S['R']
    k = S['k']
    n1 = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n1 * d1, n1 * d1)
    r1 = np.hypot(x1, y1)
    # sinc-Gaussian 模型点源（书 12-14 行；sinc 见表 1.2：sin(pi x)/(pi x)）
    with np.errstate(divide='ignore', invalid='ignore'):
        sx = np.sinc(x1 / D1)          # numpy 的 sinc 即 sin(pi x)/(pi x)
        sy = np.sinc(y1 / D1)
    pt = np.exp(-1j * k / (2 * R) * r1 ** 2) / D1 ** 2 * sx * sy \
        * np.exp(-(r1 / (4 * D1)) ** 2)
    z = np.arange(1, nplanes) * Dz / (nplanes - 1)      # 中间平面
    # 吸收边界 t（书 19-21 行：t = repmat(sg,[1 1 n])）
    sg2 = np.exp(-(x1 / (0.47 * N * d1)) ** 16) * np.exp(-(y1 / (0.47 * N * d1)) ** 16)
    t = np.repeat(sg2[:, :, None], nplanes, axis=2)
    xn, yn, Uvac = ang_spec_multi_prop(pt, wvl, d1, dn, z, t)
    # 准直：去掉残留球面波相位（书 25 行）
    Uvac = Uvac * np.exp(-1j * np.pi / (wvl * R) * (xn ** 2 + yn ** 2))
    return xn, yn, Uvac, pt, sg2, z, t


# ===== 书 Listing 9.8：湍流仿真 =====
def turbulent_simulation(S, d1, dn, N, nplanes, z, sg2, nreals=40, seed=9):
    """书 Listing 9.8：nreals 次独立湍流实现，累积互相干函数。

    返回 (xn, yn, Uout_末次, MCF2, mask, delta_scr, phz_last_all)。
    """
    wvl, Dz, D1, R, r0scrn = S['wvl'], S['Dz'], S['D1'], S['R'], S['r0scrn']
    k = S['k']
    rng = np.random.default_rng(seed)
    l0, L0 = 0.0, np.inf                # Kolmogorov PSD（书 3-4 行）
    zt = np.concatenate([[0.0], z])     # n 个平面位置
    alpha = zt / zt[-1]
    delta = (1 - alpha) * d1 + alpha * dn     # 各平面网格间距
    MCF2 = np.zeros((N, N), dtype=complex)
    n1 = np.arange(N) - N // 2
    xn0, yn0 = np.meshgrid(n1 * dn, n1 * dn)
    mask = circle_mask(xn0, yn0, S['D2'] / 2)      # 望远镜口径掩模
    Uout_last, xn, yn = None, None, None
    phz_last_all = []
    for idxreal in range(nreals):
        t = np.empty((N, N, nplanes), dtype=complex)
        for idxscr in range(nplanes):
            phz_lo, phz_hi = ft_sh_phase_screen(
                r0scrn[idxscr], N, delta[idxscr], L0, l0, rng)
            t[:, :, idxscr] = sg2 * np.exp(1j * (phz_lo + phz_hi))
            if idxscr == nplanes - 1:
                phz_last_all.append(phz_lo + phz_hi)
        xn, yn, Uout = ang_spec_multi_prop(
            _pt_source(S, d1, N), wvl, d1, dn, z, t)
        Uout = Uout * np.exp(-1j * np.pi / (wvl * R) * (xn ** 2 + yn ** 2))
        mcf, _ = corr2_ft(Uout, Uout, mask, dn)
        MCF2 = MCF2 + mcf
        Uout_last = Uout
    return xn, yn, Uout_last, MCF2, mask, delta, phz_last_all


def _pt_source(S, d1, N):
    """与 vacuum_simulation 相同的模型点源（Listing 9.7 第 12-14 行）。"""
    wvl, Dz, D1, R = S['wvl'], S['Dz'], S['D1'], S['R']
    k = S['k']
    n1 = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n1 * d1, n1 * d1)
    r1 = np.hypot(x1, y1)
    pt = np.exp(-1j * k / (2 * R) * r1 ** 2) / D1 ** 2 \
        * np.sinc(x1 / D1) * np.sinc(y1 / D1) * np.exp(-(r1 / (4 * D1)) ** 2)
    return pt


def radial_profile(Map, mean_r, rmax, nbins=120):
    """把二维量按径向分箱取平均（各向同性统计）。"""
    sel = mean_r <= rmax
    bins = np.linspace(0, rmax, nbins + 1)
    idx = np.digitize(mean_r[sel], bins) - 1
    prof = np.zeros(nbins)
    cnt = np.zeros(nbins)
    for b in range(nbins):
        m = idx == b
        if np.any(m):
            prof[b] = np.mean(Map[sel][m])
            cnt[b] = m.sum()
    centers = 0.5 * (bins[:-1] + bins[1:])
    keep = cnt > 0
    return centers[keep], prof[keep]


def main(outdir=None):
    if outdir is None:
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "..", "docs", "images", "holo", "waveprop")
    os.makedirs(outdir, exist_ok=True)

    # ================= 验证 0/1：参数与屏幕优化（Listing 9.5） =================
    S = example_setup()

    # ================= 验证 2：采样分析（Listing 9.6） =================
    SA = sampling_analysis(S)
    d1, dn, N, nplanes = SA['d1'], SA['dn'], SA['N'], S['nscr']

    # ---- 自绘 fig-9-12：采样约束图（复现书图 9.6）----
    dv1 = np.linspace(1e-4, 1.1 * S['wvl'] * S['Dz'] / SA['D2p'], 300)
    dv2 = np.linspace(1e-4, 1.1 * S['wvl'] * S['Dz'] / SA['D1p'], 300)
    D1g, DNg = np.meshgrid(dv1, dv2)
    L2N = np.log2(SA['N2'](D1g, DNg))
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    cs = ax.contourf(D1g * 1e3, DNg * 1e3, np.clip(L2N, 4, 13),
                     levels=np.arange(4, 13.5, 0.5), cmap='Greys')
    cb = fig.colorbar(cs, ax=ax); cb.set_label(r'$\log_2 N$')
    ax.plot(dv1 * 1e3, SA['dn_max'](dv1) * 1e3, 'k--', lw=2, label='约束 1')
    ax.plot(dv1 * 1e3, SA['hi3'](dv1) * 1e3, 'k-.', lw=2, label='约束 3 上界')
    ax.plot(dv1 * 1e3, np.clip(SA['lo3'](dv1), 0, None) * 1e3, 'k-.', lw=1)
    ax.plot(d1 * 1e3, dn * 1e3, 'w*', ms=14, mec='k')
    ax.set_xlim(0, 45); ax.set_ylim(0, 65)
    ax.set_xlabel(r'$\delta_1$ [mm]'); ax.set_ylabel(r'$\delta_n$ [mm]')
    ax.set_title('点源 50 km 湍流路径的采样约束（c=2, r0,sw=12.7 cm）\n'
                 r'白* = 选中 $\delta_1=\delta_n$=1 cm, $N$=512（自绘，复现书图 9.6）',
                 fontsize=10)
    ax.legend(loc='upper right', fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    fp = os.path.join(outdir, 'fig-9-12.png'); fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp); plt.close(fig)

    # ================= 验证 3：真空仿真（Listing 9.7） =================
    print('===== 验证 3：真空仿真（Listing 9.7）=====')
    xn, yn, Uvac, pt, sg2, z, t_vac = vacuum_simulation(S, d1, dn, N, nplanes)
    D2n = S['D2']
    roi = np.hypot(xn / D2n, yn / D2n) < 2.0        # |2rn/D2| < 2 的中心区
    Ivac = np.abs(Uvac) ** 2
    unc = Ivac[roi]
    print('    ROI 内辐照度：mean = %.1f W/m^2，std/mean = %.2e' % (unc.mean(), unc.std() / unc.mean()))
    ph_roi = np.angle(Uvac)[roi]
    print('    准直后 ROI 相位：RMS = %.2e rad（书图 9.7：平坦相位 + 均匀辐照度）'
          % np.sqrt(np.mean(ph_roi ** 2)))

    # ---- 自绘 fig-9-13：真空仿真结果（复现书图 9.7）----
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.0))
    ext = [2 * xn[0, 0] / D2n, 2 * xn[0, -1] / D2n,
           2 * yn[-1, 0] / D2n, 2 * yn[0, 0] / D2n]
    im0 = axs[0].imshow(Ivac, extent=ext, cmap='gray', vmin=0,
                        vmax=np.max(Ivac[np.hypot(xn / D2n, yn / D2n) < 2.5]))
    axs[0].set_title(r'(a) 辐照度 [W/m$^2$]', fontsize=10)
    plt.colorbar(im0, ax=axs[0], fraction=0.046)
    im1 = axs[1].imshow(np.angle(Uvac), extent=ext, cmap='gray', vmin=-np.pi, vmax=np.pi)
    axs[1].set_title('(b) 相位 [rad]（准直后）', fontsize=10)
    plt.colorbar(im1, ax=axs[1], fraction=0.046)
    for a in axs:
        a.set_xlabel(r'$2x_n/D_2$'); a.set_ylabel(r'$2y_n/D_2$')
    fig.suptitle('真空传播的模型点源（自绘，复现书图 9.7；N=512, '
                 r'$\delta_1=\delta_n$=1 cm, 11 平面）', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fp = os.path.join(outdir, 'fig-9-13.png'); fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp); plt.close(fig)

    # ---- 自绘 fig-9-11：相位屏方法对比（复现书图 9.3 + 典型屏=书图 9.4）----
    print('===== 验证 7：FT 屏 vs 次谐波屏（Listing 9.4 参数）=====')
    rng9 = np.random.default_rng(94)
    Ns, ds = 256, 2.0 / 256            # 2 m 屏
    r0s, L0s, l0s = 0.1, 100.0, 0.01
    nscr9 = 40
    D_ft = np.zeros((Ns, Ns)); D_sh = np.zeros((Ns, Ns))
    for _ in range(nscr9):
        ph_ft = ft_phase_screen(r0s, Ns, ds, L0s, l0s, rng9)
        plo, phi_ = ft_sh_phase_screen(r0s, Ns, ds, L0s, l0s, rng9)
        d1m, mr = str_fcn2_ft(ph_ft, ds)
        D_ft += d1m
        d2m, _ = str_fcn2_ft(plo + phi_, ds)
        D_sh += d2m
    D_ft /= nscr9; D_sh /= nscr9
    cprof, Dftp = radial_profile(D_ft, mr, 12.0 * r0s)
    _, Dshp = radial_profile(D_sh, mr, 12.0 * r0s)
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 4.2))
    plo, phi_ = ft_sh_phase_screen(r0s, Ns, ds, L0s, l0s, rng9)
    ext = [-1, 1, -1, 1]
    im = axs[0].imshow(plo + phi_, extent=ext, cmap='gray', vmin=-20, vmax=20)
    axs[0].set_title('(a) 典型次谐波相位屏 [rad]\n（D=2 m, r0=10 cm, l0=1 cm, L0=100 m）',
                     fontsize=10)
    plt.colorbar(im, ax=axs[0], fraction=0.046)
    axs[1].plot(cprof / r0s, 6.88 * (cprof / r0s) ** (5 / 3), '-', color='0.6', lw=2,
                label='理论（式 9.44）')
    axs[1].plot(cprof / r0s, Dshp, 'k-.', lw=1.5, label='Simulated SH（次谐波）')
    axs[1].plot(cprof / r0s, Dftp, 'k:', lw=1.5, label='Simulated FT（纯 FT）')
    axs[1].set_xlabel(r'$|\Delta r|/r_0$'); axs[1].set_ylabel(r'$D_\phi(|\Delta r|)$ [rad$^2$]')
    axs[1].set_ylim(0, 500); axs[1].set_xlim(0, 12)
    axs[1].set_title('(b) 40 块屏的平均结构函数（径向平均）\n（自绘，复现书图 9.3）', fontsize=10)
    axs[1].legend(fontsize=8); axs[1].grid(alpha=0.3)
    # 定量：大间隔（8-12 r0）处 SH/FT 与理论的偏差
    m8 = (cprof / r0s > 8)
    e_sh = np.max(np.abs(Dshp[m8] - 6.88 * (cprof[m8] / r0s) ** (5 / 3))) \
        / np.max(6.88 * (cprof[m8] / r0s) ** (5 / 3))
    e_ft = np.max(np.abs(Dftp[m8] - 6.88 * (cprof[m8] / r0s) ** (5 / 3))) \
        / np.max(6.88 * (cprof[m8] / r0s) ** (5 / 3))
    print('    大间隔（8-12 r0）最大相对偏差：SH = %.2f%%，FT = %.1f%%'
          % (e_sh * 100, e_ft * 100))
    print('    （低频缺失使 FT 屏大间隔结构函数远低于理论；次谐波法补救 —— 书 9.3 节）')
    fig.tight_layout()
    fp = os.path.join(outdir, 'fig-9-11.png'); fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp); plt.close(fig)

    # ================= 验证 4-6：湍流仿真 + 输出验证（Listing 9.8/9.5.5） =====
    print('===== 验证 4：湍流仿真 40 次实现（Listing 9.8）=====')
    xn, yn, Uout1, MCF2, mask, delta_scr, phz_last = turbulent_simulation(
        S, d1, dn, N, nplanes, z, sg2, nreals=40, seed=9)

    # ---- 自绘 fig-9-14：单次湍流实现（复现书图 9.8）----
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.0))
    ext = [2 * xn[0, 0] / D2n, 2 * xn[0, -1] / D2n,
           2 * yn[-1, 0] / D2n, 2 * yn[0, 0] / D2n]
    I1 = np.abs(Uout1) ** 2
    im0 = axs[0].imshow(I1, extent=ext, cmap='gray',
                        vmin=0, vmax=np.percentile(I1[np.hypot(xn / D2n, yn / D2n) < 1], 99)
                        if np.any(np.hypot(xn / D2n, yn / D2n) < 1) else np.max(I1))
    th = plt.Circle((0, 0), 1.0, color='w', fill=False, lw=1.5)
    axs[0].add_patch(th)
    axs[0].set_title(r'(a) 辐照度 [相对，白圈 = 望远镜口径]', fontsize=10)
    plt.colorbar(im0, ax=axs[0], fraction=0.046)
    im1 = axs[1].imshow(np.angle(Uout1), extent=ext, cmap='gray',
                        vmin=-np.pi, vmax=np.pi)
    axs[1].set_title('(b) 相位 [rad]（准直后）', fontsize=10)
    plt.colorbar(im1, ax=axs[1], fraction=0.046)
    for a in axs:
        a.set_xlabel(r'$2x_n/D_2$'); a.set_ylabel(r'$2y_n/D_2$')
        a.set_xlim(-2.5, 2.5); a.set_ylim(-2.5, 2.5)
    fig.suptitle('湍流传播的模型点源（40 次实现的第 1 次，自绘，复现书图 9.8）',
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fp = os.path.join(outdir, 'fig-9-14.png'); fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp); plt.close(fig)

    # ---- 验证 5：末平面屏幕的结构函数（书图 9.9）----
    print('===== 验证 5：相位屏结构函数（书图 9.9）=====')
    r0_last = S['r0scrn'][-1]
    Dacc = np.zeros((N, N))
    for phz in phz_last:
        d2m, mrN = str_fcn2_ft(phz, delta_scr[-1])
        Dacc += d2m
    Dacc /= len(phz_last)
    cprof, Dprof = radial_profile(Dacc, mrN, 12.0 * r0_last)
    err_scr = np.max(np.abs(Dprof - 6.88 * (cprof / r0_last) ** (5 / 3))) \
        / np.max(6.88 * (cprof / r0_last) ** (5 / 3))
    print('    末平面 %d 块屏（r0 = %.3f m）：与 6.88(r/r0)^(5/3) 的最大相对偏差 = %.1f%%'
          % (len(phz_last), r0_last, err_scr * 100))

    # ---- 验证 6：观测面相干因子（书图 9.10）----
    print('===== 验证 6：观测面相干因子（书图 9.10）=====')
    MCDOC2 = np.abs(MCF2) / np.abs(MCF2[N // 2, N // 2])
    n1 = np.arange(N) - N // 2
    Xr, Yr = np.meshgrid(n1 * dn, n1 * dn)
    mrdn = np.hypot(Xr, Yr)                                    # 间隔网格 [m]
    cprof, mup = radial_profile(MCDOC2, mrdn, 1.6 * S['r0sw'])
    mu_th = np.exp(-3.44 * (cprof / S['r0sw']) ** (5 / 3))
    err_mu = np.mean(np.abs(mup - mu_th))
    print('    r0,sw = %.3f cm；rn/r0 ∈ (0, 1.6] 内 |mu_sim - mu_th| 平均 = %.3f'
          % (S['r0sw'] * 100, err_mu))

    # ---- 自绘 fig-9-15：两个验证图（复现书图 9.9 + 9.10）----
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 4.2))
    axs[0].plot(cprof / r0_last, 6.88 * (cprof / r0_last) ** (5 / 3), '-',
                color='0.6', lw=2, label='理论（式 9.44）')
    axs[0].plot(cprof / r0_last, Dprof, 'k-.', lw=1.5, label='Simulated（40 块屏）')
    axs[0].set_xlabel(r'$|\Delta r|/r_0$'); axs[0].set_ylabel(r'$D_\phi(|\Delta r|)$ [rad$^2$]')
    axs[0].set_ylim(0, 500); axs[0].set_xlim(0, 12)
    axs[0].set_title('(a) 末平面屏幕平均结构函数（自绘，复现书图 9.9）', fontsize=10)
    axs[0].legend(fontsize=8); axs[0].grid(alpha=0.3)
    axs[1].plot(cprof / S['r0sw'], mu_th, '-', color='0.6', lw=2, label='theory')
    axs[1].plot(cprof / S['r0sw'], mup, 'k.-', ms=3, lw=0.8, label='simulated')
    axs[1].set_xlabel(r'$r_n/r_0$'); axs[1].set_ylabel('Coherence Factor')
    axs[1].set_ylim(0, 1); axs[1].set_xlim(0, 1.6)
    axs[1].set_title('(b) 观测面相干因子（40 次实现平均，自绘，复现书图 9.10）',
                     fontsize=10)
    axs[1].legend(fontsize=8); axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fp = os.path.join(outdir, 'fig-9-15.png'); fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp); plt.close(fig)

    print('===== 工程 1 完成 =====')
    return S, SA


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
