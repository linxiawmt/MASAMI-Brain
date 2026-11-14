# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 9 章（工程 1）=====
通过大气湍流的传播：分步束传播法 + 相位屏 + 采样 + 验证
（原书 9.1-9.6 节，Listings 9.1-9.8）

原书要点：
    分步束传播法（split-step BPM，式 9.1-9.3）：把"介质里的传播"拆成
      真空衍射算符 R 与折射（相位）算符 T——大气 n≈1，只需在第 8 章的
      真空多平面传播器（式 8.18）每个平面上乘一块相位屏。
    Kolmogorov 湍流统计（9.2.1）：惯性子区内速度/温度/折射率结构函数
      按 r^(2/3) 增长（式 9.5/9.7/9.13）；折射率 PSD（式 9.15-9.18）；
      泰勒冻结湍流假设把空间统计变成时间统计（式 9.19）。
    Rytov 微扰（9.2.2-9.2.3）：U = U0*exp(psi)，psi = chi + i*phi；
      平面波相干因子、波结构函数、Fried 参数 r0、等晕角 theta0、
      log 振幅方差（式 9.37-9.64）。平面波 Kolmogorov：
      D_phi(r) = 6.88 (r/r0)^(5/3)（式 9.44），
      r0,pw = [0.423 k^2 Cn2 Dz]^(-3/5)（式 9.42）。
    分层模型（9.2.4-9.2.5）：n 块薄相位屏匹配连续湍流的低阶矩（式 9.65），
      各屏有效相干直径 r0i（式 9.70），矩阵方程（式 9.75）+ 约束优化。
    蒙特卡洛相位屏（9.3）：FT 法（McGlamery，式 9.76-9.80）低频欠采样
      -> Lane 次谐波法补救（式 9.81，Listing 9.3）。
    采样约束（9.4）：Coy 的"光栅模糊"模型把湍流展宽折进约束 1/2：
      D1' = D1 + c*lambda*Dz/r0,rev、D2' = D2 + c*lambda*Dz/r0
      （式 9.84/9.85，c=2 约 97% 光、c=4 约 99%），约束 3 不变；
      Dz_max = min(delta1,deltan)^2*N/lambda（式 9.89）定部分传播段数。

原书勘误（本脚本按修正实现并数值验证，详见笔记）：
    (a) 式 (9.2) 印成 T = exp(-i*phi)：按时间约定 exp(-i2*pi*nu*t) 的
        抛物线方程 dU/dz = i/(2k)∇²U + i*k*δn*U，应为 exp(+i*phi)；
        书自己的 Listing 9.8（t = sg.*exp(i*phz)）用的也是 + 号。
        对零均值随机屏无统计影响，对确定性 δn（热透镜等）符号会反。
    (b) 式 (9.34) 的核函数书印 sin(kr)/(kr)（3D 场的公式）。相位是
        观测平面上的 2D 随机场，正确的逆变换核是 J1：
        Phi(kappa) = 1/(4*pi*kappa) * ∫ J1(kappa*r) D'(r) r dr，
        代入式 (9.44) 恰得 0.49（式 9.49）；sin 核（加阻尼数值积分）
        给 ~0.35（习题 3 实算，见 ch9_problems.py）。
    (c) 式 (9.48) 系数书印 2*pi^2*k^2*Dz，应为 2*pi*k^2*Dz：
        与 (9.42)+(9.16) 联立应给 (9.49) 的 0.49，2*pi*0.033/0.423 = 0.490
        （对）；2*pi^2 给 1.545（错，差 pi 倍）。
    (d) 正文 9.5.4 说"重复 40 次"，Listing 9.8 却印 nreals = 20（两处
        不一致）；本脚本取 20（与 Listing 一致）。
    (e) Listing 9.8 第 31 行续行符书印 4 个点 "...."，应为 "..."。
    (f) 式 (9.11)/(9.12)（dn 的微分式）书印丢因子：按 n-1 = 7.99e-5*P/T
        微分应为 dn = 7.99e-5*(dP/T - P*dT/T^2)；dP≈0 时
        dn = -7.99e-5*P*dT/T^2，与 dθ 成正比（差一个缓变因子与整体符号，
        都被吸收进正定的结构参数 Cn^2），式 (9.13) 之后不受影响。

实现说明（与书代码的对应关系）：
    - ft2/ift2：复用 fourier/common/base.py 的 ft2_centered 与
      waveprop/ch6 的 ift2（书 Listing 2.5/2.6）。
    - str_fcn2_ft / corr2_ft：书第 3 章 Listing 3.7/3.10（ch3 PDF 缺失），
      按 9.5.5 节用法重写：str_fcn2_ft 按定义算无偏"逐对差分"结构函数
      D(rho) = <(phi(r+rho)-phi(r))^2>（填零到 2N 的 FFT 线性相关 +
      重叠样本数归一）。不能用 2(var-cov) 的"方差-协方差"捷径：次谐波
      屏含波长超过屏幕的大幅低频模，全局方差不等于逐对平均的均方。
      corr2_ft 用相关定理算口径内互相干函数（口径 << 网格，循环相关
      即线性相关），第二个返回值是各像素的间隔模长（供径向分箱）。
    - ft_phase_screen 合成只取 IFT 实部：cn = (randn+i*randn)*sqrt(PSD)*df
      的 |cn|^2 期望是式 (9.79) 的 2 倍，Re[.] 恰好带来 1/2，两者相消
      ——Re/Im 各自都是方差正确的独立相位屏（书 Listing 9.2 第 20 行）。
    - 随机数用 numpy default_rng(seed) 播种，结果可复现（书未播种）。
    - 屏 r0 的约束优化：书用 fmincon（Optimization Toolbox）；本脚本是
      2 x nscr 矩阵的有界线性最小二乘，用 scipy.optimize.minimize
      （L-BFGS-B，与 fmincon 同为"目标函数+边界"的一般框架），
      初值/边界与书一致，回代校验 r0,sw 与 Rytov 数。

用法（命令行分步跑，单步 < 60 s；输出 out_*.png + out_*.npz 到 outdir）：
    python ch9_turbulence.py screens [outdir]   # 工程 A：FT 屏 vs 次谐波屏
                                                # 结构函数（书图 9.3 复现）
    python ch9_turbulence.py main [outdir]      # 书例（9.5 节）：50 km、
                                                # Cn2=1e-16、11 屏、20 次实现
    python ch9_turbulence.py p7 [outdir]        # 习题 7：100 km、Cn2=1e-17、
                                                # 7 屏、20 次实现（D2=0.5 m 示例规模）
    python ch9_turbulence.py weak [outdir]      # 弱湍流对照：50 km、Cn2=1e-17
    python ch9_turbulence.py trend [outdir]     # 汇总三例：闪烁随 Rytov 方差变化
依赖：numpy、scipy、matplotlib；base（fourier/common）、ch6_fresnel_prop。
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'dejavusans'   # 数学字体用 DejaVu（含 U+2212）

# ---- 复用第 2/6 章的基本算子（同仓库） ----
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'ch6'))
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'fourier', 'common'))
from base import ft2_centered as ft2                       # noqa: E402  书 Listing 2.5
from ch6_fresnel_prop import ift2                          # noqa: E402  书 Listing 2.6


# =====================================================================
# Listing 9.1：ang_spec_multi_prop —— 分步束传播法（式 9.3）
# =====================================================================
def ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t):
    """角谱多平面传播（书 Listing 9.1）。

    t: (N, N, n) 相互作用数组——吸收边界（真空）或 吸收边界*exp(i*相位屏)
    （湍流）；t[:,:,0] 作用于源面，t[:,:,i+1] 作用于第 i 次部分传播后。
    T=1（全 1 数组）时即第 8 章的真空多平面传播。
    注意：书实现在观测面也乘吸收边界（ROI 内 sg≈1，无实际影响）。"""
    N = Uin.shape[0]
    n1 = np.arange(N) - N // 2                    # MATLAB (-N/2 : N/2-1)
    k = 2 * np.pi / wvl
    # 超高斯吸收边界（第 8 章）：sg = exp(-nsq^8/w^16) = exp(-(nsq/w^2)^8)
    nsq = n1[None, :] ** 2 + n1[:, None] ** 2
    w = 0.47 * N
    sg = np.exp(-(nsq / w ** 2) ** 8)

    z = np.concatenate(([0.0], np.asarray(z, dtype=float)))
    n = len(z)
    Delta_z = np.diff(z)
    alpha = z / z[-1]                             # 式 (8.8)：网格间距线性插值
    delta = (1 - alpha) * delta1 + alpha * deltan
    m = delta[1:] / delta[:-1]
    x1 = n1[None, :] * delta[0]
    y1 = n1[:, None] * delta[0]
    r1sq = x1 ** 2 + y1 ** 2
    Q1 = np.exp(1j * k / 2 * (1 - m[0]) / Delta_z[0] * r1sq)
    Uin = Uin * Q1 * t[:, :, 0]
    for idx in range(n - 1):
        deltaf = 1 / (N * delta[idx])
        fX = n1[None, :] * deltaf
        fY = n1[:, None] * deltaf
        fsq = fX ** 2 + fY ** 2
        Z = Delta_z[idx]
        Q2 = np.exp(-1j * np.pi ** 2 * 2 * Z / m[idx] / k * fsq)
        Uin = sg * t[:, :, idx + 1] * ift2(Q2 * ft2(Uin / m[idx], delta[idx]), deltaf)
    xn = n1[None, :] * delta[-1]
    yn = n1[:, None] * delta[-1]
    rnsq = xn ** 2 + yn ** 2
    Q3 = np.exp(1j * k / 2 * (m[-1] - 1) / (m[-1] * Z) * rnsq)   # Z = 末段距离
    Uout = Q3 * Uin
    return xn, yn, Uout


# =====================================================================
# Listing 9.2 / 9.3：相位屏生成（FT 法 + 次谐波）
# =====================================================================
def _mvk_phase_psd(f, r0, L0, l0):
    """修正 von Karman 相位 PSD（式 9.51 的 f 版，式 9.52：0.023 r0^-5/3）。
    l0=0 / L0=inf 时分别退化为无内尺度 / Kolmogorov 谱。"""
    fm = 5.92 / l0 / (2 * np.pi) if l0 > 0 else np.inf   # kappa_m = 5.92/l0
    f0 = 1.0 / L0 if np.isfinite(L0) else 0.0            # kappa_0 = 2*pi/L0
    with np.errstate(divide='ignore', invalid='ignore'):
        psd = 0.023 * r0 ** (-5.0 / 3) * np.exp(-(f / fm) ** 2) \
            / (f ** 2 + f0 ** 2) ** (11.0 / 6)
    return psd


def ft_phase_screen(r0, N, delta, L0, l0, rng):
    """FT 法相位屏（书 Listing 9.2）。
    E{|cn|^2} = PSD*df^2（式 9.79/9.80）；(randn+i*randn) 的模方期望
    是 2 倍、取实部又带来 1/2，恰好相消——real(IFT) 即方差正确的屏。"""
    del_f = 1 / (N * delta)
    fx = (np.arange(N) - N // 2) * del_f
    FX, FY = np.meshgrid(fx, fx)          # 与 MATLAB meshgrid(fx) 一致
    f = np.hypot(FX, FY)
    psd = _mvk_phase_psd(f, r0, L0, l0)
    psd[N // 2, N // 2] = 0               # 零频（piston）置零
    cn = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) \
        * np.sqrt(psd) * del_f
    return np.real(ift2(cn, 1.0))


def ft_sh_phase_screen(r0, N, delta, L0, l0, rng):
    """FT 法 + 次谐波（书 Listing 9.3，Lane et al.）。
    低频补偿屏：p=1,2,3 三层网格（间距 1/(3^p*L)）各取 3×3 个最低频
    分量，按 PSD 抽系数直接叠加平面波（频率少，不必 FFT）。
    返回 (phz_lo, phz_hi)，总屏 = phz_lo + phz_hi。"""
    D = N * delta
    phz_hi = ft_phase_screen(r0, N, delta, L0, l0, rng)
    x = (np.arange(N) - N // 2) * delta
    X, Y = np.meshgrid(x, x)
    phz_lo = np.zeros((N, N), dtype=complex)
    for p in range(1, 4):
        del_f = 1 / (3 ** p * D)
        f1 = np.array([-1.0, 0.0, 1.0]) * del_f
        FX, FY = np.meshgrid(f1, f1)      # 3×3 频点（MATLAB 线性索引同款）
        f = np.hypot(FX, FY)
        psd = _mvk_phase_psd(f, r0, L0, l0)
        psd[1, 1] = 0                     # (0,0) 频点置零（书 PSD_phi(2,2)）
        cn = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) \
            * np.sqrt(psd) * del_f
        SH = np.zeros((N, N), dtype=complex)
        for ii in range(3):
            for jj in range(3):
                SH += cn[ii, jj] * np.exp(1j * 2 * np.pi
                                          * (FX[ii, jj] * X + FY[ii, jj] * Y))
        phz_lo += SH
    phz_lo = np.real(phz_lo) - np.real(phz_lo).mean()
    return phz_lo, phz_hi


# =====================================================================
# Listing 3.7 / 3.10 的重写：结构函数与互相干函数
# =====================================================================
def str_fcn2_ft(phz, delta):
    """无偏二维相位结构函数（书 Listing 3.7 重写，ch3 PDF 缺失）。
    D(rho) = <(phi(r+rho)-phi(r))^2>，填零到 2N 的 FFT 线性相关 +
    重叠样本数归一。返回 (D, rgrid)：D 是 rho_x,rho_y>=0 象限、
    零间隔在 [0,0]，D[0, :] 即沿 x 轴的径向切片（r = rgrid[0, :]）。
    不用 2(var-cov) 捷径——次谐波屏低频模波长超过屏幕，全局方差
    不等于逐对平均均方（小间隔高估数倍）。"""
    N = phz.shape[0]
    N2 = 2 * N
    P0 = np.zeros((N2, N2)); P0[:N, :N] = phz
    W0 = np.zeros((N2, N2)); W0[:N, :N] = 1.0
    Q0 = np.zeros((N2, N2)); Q0[:N, :N] = phz ** 2
    FP = np.fft.rfft2(P0)
    FW = np.fft.rfft2(W0)
    FQ = np.fft.rfft2(Q0)
    C = np.fft.irfft2(np.conj(FP) * FP, s=(N2, N2))      # sum phi(j)phi(j+m)
    S1 = np.fft.irfft2(np.conj(FQ) * FW, s=(N2, N2))     # sum phi(j)^2 W(j+m)
    S2 = np.fft.irfft2(FQ * np.conj(FW), s=(N2, N2))     # sum W(j) phi(j+m)^2
    Np = np.fft.irfft2(np.conj(FW) * FW, s=(N2, N2))     # 重叠像素数
    Dfull = (S1 + S2 - 2 * C) / np.maximum(Np, 1)
    # wrap 序里正移位就在 [0, N)：直接取左上象限（与 MATLAB 版一致）
    Dq = Dfull[:N, :N]
    ii = np.arange(N)[:, None] * delta
    jj = np.arange(N)[None, :] * delta
    return Dq, np.hypot(ii, jj)


def corr2_ft(U1, U2, mask, delta):
    """口径内互相干函数（书 Listing 3.10 重写）：
    MCF(rho) = sum_r U1(r) U2*(r+rho) mask(r) mask(r+rho)。
    口径直径 << 网格边长时循环相关 == 线性相关（无卷绕污染）。
    返回 fftshift 后的居中数组（零间隔在 [N//2, N//2]，与书一致）+
    各像素间隔模长 rgrid。"""
    N = U1.shape[0]
    A = np.fft.fft2(U1 * mask)
    B = np.fft.fft2(U2 * mask)      # 共轭在频域做：corr 定理 IFT[conj(A)B]
    mcf = np.fft.fftshift(np.fft.ifft2(np.conj(A) * B))
    n1 = np.arange(N) - N // 2
    rr = n1[:, None] * delta
    cc = n1[None, :] * delta
    return mcf, np.hypot(rr, cc)


# =====================================================================
# 9.2.5 / Listing 9.5：分层大气矩匹配 + 屏 r0 约束优化
# =====================================================================
def solve_screen_r0s(r0sw_val, rytov_val, nscr, k, Dz, rmax=0.1):
    """由目标 (r0,sw, sigma_chi^2,sw) 反解 nscr 块屏的 r0i（式 9.72/9.74）。
    与书 Listing 9.5 的 fmincon 同框架：目标 ||A X - b||^2（X = r0i^-5/3），
    下界 0（r0->inf = 无湍流），上界使每屏 Rytov 贡献 <= rmax
    （Martin & Flatté 指南）；权重为 0 的屏（路径两端）上界放宽到
    r0 = 50 m（书 x2(A(2,:)==0) = 50^(-5/3)）。"""
    alpha = np.arange(nscr) / (nscr - 1)
    A = np.vstack([alpha ** (5.0 / 3),
                   (1 - alpha) ** (5.0 / 6) * alpha ** (5.0 / 6)])
    b = np.array([r0sw_val ** (-5.0 / 3),
                  rytov_val / 1.33 * (k / Dz) ** (5.0 / 6)])
    with np.errstate(divide='ignore', invalid='ignore'):
        ub = rmax / 1.33 * (k / Dz) ** (5.0 / 6) / A[1]
    ub[A[1] == 0] = 50 ** (-5.0 / 3)
    x0 = np.full(nscr, (nscr / 3 * r0sw_val) ** (-5.0 / 3))
    fun = lambda X: float(np.sum((A @ X - b) ** 2))
    res = minimize(fun, x0, bounds=list(zip(np.zeros(nscr), ub)),
                   method='L-BFGS-B',
                   options={'maxiter': 2000, 'ftol': 1e-18, 'gtol': 1e-14})
    X = res.x
    with np.errstate(divide='ignore', invalid='ignore'):
        r0scrn = X ** (-3.0 / 5)
    r0scrn[~np.isfinite(r0scrn)] = 1e6
    bp = A @ X
    chk = (bp[0] ** (-3.0 / 5), bp[1] * 1.33 * (Dz / k) ** (5.0 / 6))
    return r0scrn, alpha, chk


def atmos_setup(D2, wvl, Dz, Cn2, nscr, rmax=0.1):
    """书例（Listing 9.5）的几何 + 大气参数 + 屏 r0 优化。"""
    k = 2 * np.pi / wvl
    DROI = 4 * D2                       # 观测面感兴趣区直径（均匀照亮）
    D1 = wvl * Dz / DROI                # sinc 点源模型中央瓣全宽 2*D1
    R = Dz                              # 波前曲率半径
    # 常数 Cn2 的解析积分：式 (9.43)/(9.42)/(9.64)
    r0sw = (0.423 * k ** 2 * Cn2 * 3 / 8 * Dz) ** (-3.0 / 5)
    r0pw = (0.423 * k ** 2 * Cn2 * Dz) ** (-3.0 / 5)
    p = np.linspace(0, Dz, 2001)
    dp = p[1] - p[0]
    rytov = 0.563 * k ** (7.0 / 6) * np.sum(Cn2 * (1 - p / Dz) ** (5.0 / 6)
                                            * p ** (5.0 / 6)) * dp
    r0scrn, alpha, chk = solve_screen_r0s(r0sw, rytov, nscr, k, Dz, rmax)
    return dict(D2=D2, wvl=wvl, k=k, Dz=Dz, Cn2=Cn2, DROI=DROI, D1=D1, R=R,
                r0sw=r0sw, r0pw=r0pw, rytov=rytov, nscr=nscr,
                r0scrn=r0scrn, alpha=alpha, chk=chk)


# =====================================================================
# 9.5.3-9.5.5：真空/湍流仿真 campaign + 验证统计
# =====================================================================
def _point_source(N, delta1, wvl, k, D1, R):
    """sinc-Gaussian 点源模型（Listing 9.7 的 pt）。"""
    n1 = np.arange(N) - N // 2
    x1 = n1[None, :] * delta1
    y1 = n1[:, None] * delta1
    r1 = np.hypot(x1, y1)
    return (np.exp(-1j * k / (2 * R) * r1 ** 2) / D1 ** 2
            * np.sinc(x1 / D1) * np.sinc(y1 / D1)
            * np.exp(-(r1 / (4 * D1)) ** 2))


def _absorbing(x1g, y1g, N, d1):
    """超高斯吸收边界（Listing 9.7 的 sg，分离形式）。"""
    return (np.exp(-(x1g / (0.47 * N * d1)) ** 16)
            * np.exp(-(y1g / (0.47 * N * d1)) ** 16))


def remove_tilt(U, X, Y, mask, span=64.0, n_round=5):
    """短曝光近似：估计并去掉每块实现的总倾斜（Zernike tilt 模）。
    在口径内找 (bx, by) 使相干和 |sum U e^{-i(bx x + by y)}| 最大
    （= 口径上最优拟合平面），分级粗扫细化，返回（去倾斜场, 倾斜）。"""
    idx = np.where(mask > 0)
    Xs, Ys, Us = X[idx], Y[idx], U[idx]

    def g(bx, by):
        return abs(np.sum(Us * np.exp(-1j * (bx * Xs + by * Ys))))

    bx = by = 0.0
    step = span
    for _ in range(n_round):                # 分级搜索：64 -> 64/8^4 rad/m
        best = g(bx, by)
        for dbx in np.arange(-step, step + 1e-9, step / 4):
            for dby in np.arange(-step, step + 1e-9, step / 4):
                val = g(bx + dbx, by + dby)
                if val > best:
                    best, bx, by = val, bx + dbx, by + dby
        step /= 8.0
    a = np.angle(np.sum(Us * np.exp(-1j * (bx * Xs + by * Ys))))
    return U * np.exp(-1j * (a + bx * X + by * Y)), (bx, by)


def run_campaign(tag, outdir, D2=0.5, Dz=50e3, Cn2=1e-16, nscr=11,
                 nreals=20, seed=9, l0=0.0, L0=np.inf):
    """完整 campaign：建立几何/大气 -> 采样 -> 真空 -> nreals 次湍流实现
    -> 相干因子（长曝光 = 原始；短曝光近似 = 去倾斜）+ 闪烁统计
    + 中程屏结构函数验证。返回结果字典并存 npz。"""
    os.makedirs(outdir, exist_ok=True)
    st = atmos_setup(D2=D2, wvl=1e-6, Dz=Dz, Cn2=Cn2, nscr=nscr)
    wvl, k, D1, R = st['wvl'], st['k'], st['D1'], st['R']
    # ---- 采样分析（式 9.86-9.90，c=2）----
    c = 2
    D1p = D1 + c * wvl * Dz / st['r0sw']            # r0,rev = r0,sw（常数 Cn2）
    D2p = D2 + c * wvl * Dz / st['r0sw']
    delta1 = 2 * D1 / 5                              # 点源中央瓣 5 点采样
    deltan = min(D2 / 50, (wvl * Dz - D2p * delta1) / D1p)   # 约束 1 + 口径 50 点
    Nmin = (D1p / (2 * delta1) + D2p / (2 * deltan)
            + wvl * Dz / (2 * delta1 * deltan))      # 约束 2（式 9.87）
    N = int(2 ** np.ceil(np.log2(Nmin)))
    zmax = min(delta1, deltan) ** 2 * N / wvl        # 式 (9.89)
    nmin_planes = int(np.ceil(Dz / zmax)) + 1        # 式 (9.90)
    print('===== %s：几何/大气/采样 =====' % tag)
    print('    Dz=%.0f km, Cn2=%.1e, D2=%.2f m, nscr=%d' % (Dz / 1e3, Cn2, D2, nscr))
    print('    r0,pw = %.4f m, r0,sw = %.4f m, sigma_chi^2,sw = %.4f（式 9.42/9.43/9.64）'
          % (st['r0pw'], st['r0sw'], st['rytov']))
    print('    约束 1: deltan <= %.4f m; 取 delta1=%.4f m, deltan=%.4f m'
          % ((wvl * Dz - D2p * delta1) / D1p, delta1, deltan))
    print('    约束 2: N >= %.1f -> N = 2^%d = %d; Dz_max = %.1f km -> 最少 %d 平面'
          % (Nmin, np.ceil(np.log2(Nmin)), N, zmax / 1e3, nmin_planes))
    print('    屏 r0 优化回代：r0,sw = %.4f m（目标 %.4f）、sigma^2 = %.4f（目标 %.4f）'
          % (st['chk'][0], st['r0sw'], st['chk'][1], st['rytov']))

    # ---- 网格与源 ----
    n = nscr                                         # 平面数（每屏一个平面）
    z = (np.arange(1, n)) * Dz / (n - 1)             # n-1 次部分传播
    x1 = (np.arange(N) - N // 2) * delta1
    X1g, Y1g = np.meshgrid(x1, x1)
    pt = _point_source(N, delta1, wvl, k, D1, R)
    sg2 = _absorbing(X1g, Y1g, N, delta1)
    t_vac = np.repeat(sg2[:, :, None], n, axis=2)
    # 真空仿真（Listing 9.7）
    xn, yn, Uvac = ang_spec_multi_prop(pt, wvl, delta1, deltan, z, t_vac)
    Uvac = Uvac * np.exp(-1j * np.pi / (wvl * R) * (xn ** 2 + yn ** 2))   # 准直
    Ivac = np.abs(Uvac) ** 2
    # 均匀性度量取平台内部（|x|<=0.75*DROI/2）：sinc 傅里叶卷积把 rect 边缘
    # 平滑掉 ~0.1 m，平台边缘不属于"均匀照亮区"（书图 9.7a 同样只在
    # 中部平坦）
    roi = ((np.abs(xn) <= 0.75 * st['DROI'] / 2)
           & (np.abs(yn) <= 0.75 * st['DROI'] / 2))
    phm = roi & (Ivac > 0.05 * Ivac.max())
    ph_plateau = np.angle(Uvac[phm])
    ph_plateau = ph_plateau - ph_plateau.mean()      # 去掉全局 piston（含 1/i 的 -pi/2）
    print('    真空平台（|x|<=0.75*DROI/2）：I mean = %.1f W/m^2, std/mean = %.2e；'
          '准直后相位 std = %.2e rad（去 piston，应近 0）'
          % (Ivac[roi].mean(), Ivac[roi].std() / Ivac[roi].mean(),
             ph_plateau.std()))

    # ---- 湍流实现（Listing 9.8）----
    zt = np.concatenate(([0.0], z))
    alpha = zt / zt[-1]
    delta = (1 - alpha) * delta1 + alpha * deltan    # 各平面网格间距
    mask = ((xn / D2) ** 2 + (yn / D2) ** 2 <= 1.0).astype(float)   # 望远镜口径
    Xn2d = np.broadcast_to(xn, (N, N))               # 去倾斜用的满二维坐标
    Yn2d = np.broadcast_to(yn, (N, N))
    rng = np.random.default_rng(seed)
    MCF2 = np.zeros((N, N), dtype=complex)
    MCF2t = np.zeros((N, N), dtype=complex)
    s_chi = np.zeros((N, N)); s_chi2 = np.zeros((N, N))          # log 振幅
    s_I = np.zeros((N, N)); s_I2 = np.zeros((N, N))              # 辐照度
    Dscr_mid = np.zeros(N); Dscr_last = np.zeros(N)              # 沿 x 轴切片
    mid = n // 2
    r0_mid = st['r0scrn'][mid]
    I_t0 = ph_t0 = None
    for ireal in range(nreals):
        phz = np.zeros((N, N, n))
        for iscr in range(n):
            lo, hi = ft_sh_phase_screen(st['r0scrn'][iscr], N, delta[iscr],
                                        L0, l0, rng)
            phz[:, :, iscr] = lo + hi
        t_turb = t_vac * np.exp(1j * phz)            # 勘误 a：exp(+i*phi)
        _, _, Uout = ang_spec_multi_prop(pt, wvl, delta1, deltan, z, t_turb)
        Uout = Uout * np.exp(-1j * np.pi / (wvl * R) * (xn ** 2 + yn ** 2))
        mcf, _ = corr2_ft(Uout, Uout, mask, deltan)
        MCF2 += mcf
        Ut, _ = remove_tilt(Uout, Xn2d, Yn2d, mask)
        mcf2, _ = corr2_ft(Ut, Ut, mask, deltan)
        MCF2t += mcf2
        chi = np.log(np.abs(Uout) / np.abs(Uvac))    # log 振幅扰动（式 9.62）
        It = np.abs(Uout) ** 2
        s_chi += chi; s_chi2 += chi ** 2
        s_I += It; s_I2 += It ** 2
        Dq, _ = str_fcn2_ft(phz[:, :, mid], delta[mid])
        Dscr_mid += Dq[0, :]
        Dq, _ = str_fcn2_ft(phz[:, :, -1], delta[-1])
        Dscr_last += Dq[0, :]
        if ireal == 0:
            I_t0, ph_t0 = It, np.angle(Uout)
        if (ireal + 1) % 5 == 0:
            print('    实现 %d/%d 完成' % (ireal + 1, nreals))

    # ---- 验证统计 ----
    MCDOC2 = np.abs(MCF2) / np.real(MCF2[N // 2, N // 2])        # 长曝光（含倾斜）
    MCDOC2t = np.abs(MCF2t) / np.real(MCF2t[N // 2, N // 2])     # 短曝光近似
    rr = np.hypot(xn, yn)
    rb = np.linspace(0.02, 1.5, 60) * st['r0sw']
    prof, proft = [], []
    for r1_ in rb:
        selm = (rr >= r1_ - deltan / 2) & (rr < r1_ + deltan / 2)
        prof.append(np.mean(MCDOC2[selm]))
        proft.append(np.mean(MCDOC2t[selm]))
    prof, proft = np.array(prof), np.array(proft)
    muth = np.exp(-3.44 * (rb / st['r0sw']) ** (5.0 / 3))    # 式 (9.32)+(9.44)
    dev = np.mean(np.abs(prof - muth))
    sel_fit = (rb >= 0.3 * st['r0sw']) & (rb <= 1.0 * st['r0sw'])
    x_fit = (rb[sel_fit] / st['r0sw']) ** (5.0 / 3)
    cf = np.polyfit(x_fit, np.log(np.maximum(prof[sel_fit], 1e-12)), 1)
    r0_meas = st['r0sw'] * (3.44 / (-cf[0])) ** (3.0 / 5)    # 由斜率反推等效 r0
    # 闪烁：log 振幅方差（式 9.62）与闪烁指数
    var_chi = (s_chi2 / nreals - (s_chi / nreals) ** 2)
    sigI = s_I2 / nreals / (s_I / nreals) ** 2 - 1
    # 中程屏结构函数（书图 9.9 同类验证）：沿 x 轴切片 + 2-8 r0 段斜率
    r_ax = np.arange(N) * delta[mid]
    Dmid = Dscr_mid / nreals
    Dlast = Dscr_last / nreals
    sl_semid = (r_ax >= 2 * r0_mid) & (r_ax <= 8 * r0_mid)
    cf2 = np.polyfit(np.log(r_ax[sl_semid]), np.log(np.maximum(Dmid[sl_semid], 1e-9)), 1)
    sl_near = (r_ax >= 0.5 * r0_mid) & (r_ax <= 2 * r0_mid)
    cf2n = np.polyfit(np.log(r_ax[sl_near]), np.log(np.maximum(Dmid[sl_near], 1e-9)), 1)
    D_r0 = np.interp(r0_mid, r_ax, Dmid)
    print('    相干因子（长曝光）：0<r/r0<=1.5 内 |mu_sim-mu_th| 平均 = %.3f；'
          'ln(mu)~(r/r0)^(5/3) 拟合斜率 = %.2f（理论 -3.44），等效 r0 = %.4f m'
          % (dev, cf[0], r0_meas))
    print('    去倾斜（短曝光近似）：r=r0,sw 处 mu = %.3f（长曝光 %.3f，'
          '理论 e^{-3.44} = 0.032）'
          % (proft[np.argmin(np.abs(rb - st['r0sw']))],
             prof[np.argmin(np.abs(rb - st['r0sw']))]))
    print('    闪烁：口径内 <sigma_chi^2> = %.4f（预测 %.4f，式 9.64）；'
          '<sigma_I^2> = %.4f（弱起伏理论 4*sigma_chi^2 = %.4f）'
          % (var_chi[mask > 0].mean(), st['rytov'], sigI[mask > 0].mean(),
             4 * st['rytov']))
    print('    中程屏（r0 = %.3f m）结构函数：0.5-2 r0 段斜率 = %.3f、'
          '2-8 r0 段斜率 = %.3f（理论 5/3 = 1.667）；D(r0) = %.2f（理论 6.88）'
          % (r0_mid, cf2n[0], cf2[0], D_r0))

    res = dict(tag=tag, st=st, N=N, delta1=delta1, deltan=deltan, n=n, z=z,
               nreals=nreals, rb=rb, prof=prof, proft=proft, dev=dev,
               r0_meas=r0_meas, var_chi=float(var_chi[mask > 0].mean()),
               sigI=float(sigI[mask > 0].mean()), rytov=st['rytov'],
               Ivac=Ivac, Uvac_ph=np.angle(Uvac), xn=xn, yn=yn,
               I_t0=I_t0, ph_t0=ph_t0, Dmid=Dmid, Dlast=Dlast,
               r_ax=r_ax, r0_mid=r0_mid, slope_mid=cf2[0], D_r0=D_r0)
    np.savez_compressed(os.path.join(outdir, 'out_%s.npz' % tag),
                        rb=rb, prof=prof, proft=proft, dev=dev,
                        r0_meas=r0_meas, var_chi=var_chi[mask > 0].mean(),
                        sigI=sigI[mask > 0].mean(), rytov=st['rytov'],
                        r0sw=st['r0sw'], r0pw=st['r0pw'],
                        r0scrn=st['r0scrn'], alpha=st['alpha'],
                        Dmid=Dmid, Dlast=Dlast, r_ax=r_ax, r0_mid=r0_mid,
                        slope_mid=cf2[0], D_r0=D_r0)
    # 图：真空/湍流场（一次代表性实现）+ 分层屏 + 相干因子
    _plot_fields(res, outdir)
    _plot_layering(res, outdir)
    _plot_coh(res, outdir)
    return res


def _plot_fields(res, outdir):
    """真空 vs 湍流的辐照度/相位（对照书图 9.7/9.8）。"""
    st = res['st']
    xn, yn = res['xn'], res['yn']
    ext = (xn[0, 0] / st['D2'], xn[0, -1] / st['D2'],
           yn[0, 0] / st['D2'], yn[-1, 0] / st['D2'])
    fig, axs = plt.subplots(2, 2, figsize=(10.8, 8.8))
    im = axs[0, 0].imshow(res['Ivac'], extent=ext, origin='lower', cmap='gray')
    fig.colorbar(im, ax=axs[0, 0], label=r'Irradiance [W/m$^2$]')
    axs[0, 0].set_title('(a) 真空辐照度（书图 9.7a）', fontsize=11)
    im = axs[0, 1].imshow(res['Uvac_ph'], extent=ext, origin='lower', cmap='gray')
    fig.colorbar(im, ax=axs[0, 1], label='Phase [rad]')
    axs[0, 1].set_title('(b) 真空相位（准直后，书图 9.7b）', fontsize=11)
    im = axs[1, 0].imshow(res['I_t0'], extent=ext, origin='lower', cmap='gray')
    fig.colorbar(im, ax=axs[1, 0], label=r'Irradiance [W/m$^2$]')
    th = np.linspace(0, 2 * np.pi, 256)
    axs[1, 0].plot(np.cos(th), np.sin(th), 'w-', lw=1.2)
    axs[1, 0].set_xlim(-2, 2); axs[1, 0].set_ylim(-2, 2)
    axs[1, 0].set_title('(c) 湍流辐照度（书图 9.8a，白圈=口径）', fontsize=11)
    im = axs[1, 1].imshow(res['ph_t0'], extent=ext, origin='lower', cmap='gray')
    fig.colorbar(im, ax=axs[1, 1], label='Phase [rad]')
    axs[1, 1].plot(np.cos(th), np.sin(th), 'w-', lw=1.2)
    axs[1, 1].set_xlim(-2, 2); axs[1, 1].set_ylim(-2, 2)
    axs[1, 1].set_title('(d) 湍流相位（准直后，书图 9.8b）', fontsize=11)
    for ax in axs.ravel():
        ax.set_xlabel(r'$2x_n/D_2$'); ax.set_ylabel(r'$2y_n/D_2$')
    fig.suptitle('图 9-13（自绘仿真）点源 %s 真空/湍流传播（书 9.5.3-9.5.4 节，'
                 'N=512, 11 平面, 20 次实现）' % ('50 km'), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(outdir, 'out_vac_turb_fields.png')
    fig.savefig(fp, dpi=150); plt.close(fig)
    print('    已保存 %s' % fp)


def _plot_layering(res, outdir):
    """屏 r0 分布与矩匹配权重（书 9.2.5 节 + Listing 9.5 结果）。"""
    st = res['st']
    al = st['alpha']
    r0s = st['r0scrn'].copy()
    fig, axs = plt.subplots(1, 2, figsize=(11.4, 4.3))
    axs[0].stem(al, r0s, basefmt=' ')
    axs[0].set_yscale('log')
    axs[0].set_ylim(0.05, 200)
    axs[0].set_xlabel(r'$\alpha_i = z_i/\Delta z$')
    axs[0].set_ylabel(r'$r_{0i}$ [m]')
    axs[0].set_title('(a) 约束优化解出的 11 块屏强度（Listing 9.5）\n'
                     '（两端权重为 0 的屏被上界/初值钉住）', fontsize=11)
    axs[0].grid(alpha=0.3)
    a2 = np.linspace(0, 1, 200)
    axs[1].plot(a2, a2 ** (5.0 / 3), 'k-',
                label=r'$\alpha^{5/3}$（$r_{0,sw}$ 行，式 9.72）')
    axs[1].plot(a2, (a2 * (1 - a2)) ** (5.0 / 6), 'k--',
                label=r'$[\alpha(1-\alpha)]^{5/6}$（Rytov 行，式 9.74）')
    axs[1].plot(a2, (1 - a2) ** (5.0 / 6), 'k-.',
                label=r'$(1-\alpha)^{5/6}$（$\sigma^2_{\chi,pw}$ 行，式 9.73）')
    axs[1].set_xlabel(r'$\alpha_i$'); axs[1].legend(fontsize=8)
    axs[1].set_title('(b) 矩匹配方程的权重函数', fontsize=11)
    axs[1].grid(alpha=0.3)
    fig.suptitle('图 9-12（自绘仿真）分层大气：屏 r0 解与矩匹配权重', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fp = os.path.join(outdir, 'out_layered_screens.png')
    fig.savefig(fp, dpi=150); plt.close(fig)
    print('    已保存 %s' % fp)


def _plot_coh(res, outdir):
    """相干因子：理论 / 长曝光 / 去倾斜短曝光近似 + 末屏结构函数。"""
    st = res['st']
    fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.8))
    ax = axs[0]
    ax.plot(res['rb'] / st['r0sw'],
            np.exp(-3.44 * (res['rb'] / st['r0sw']) ** (5.0 / 3)),
            'k-', lw=2,
            label=r'理论 $\exp[-3.44(r/r_{0,sw})^{5/3}]$（式 9.32+9.44）')
    ax.plot(res['rb'] / st['r0sw'],
            np.exp(res['rytov'] - 3.44 * (res['rb'] / st['r0sw']) ** (5.0 / 3)),
            color='gray', ls=':', lw=1.2,
            label=r'含闪烁幅度修正 $\mathrm{e}^{\sigma_\chi^2}$')
    ax.plot(res['rb'] / st['r0sw'], res['prof'], 'k.', ms=6,
            label='仿真（长曝光，含倾斜）')
    ax.plot(res['rb'] / st['r0sw'], res['proft'], 'b^', ms=5, mfc='none',
            label='仿真（去倾斜 = 短曝光近似）')
    ax.set_xlabel(r'$r_n / r_{0,sw}$'); ax.set_ylabel('相干因子')
    ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title('(a) 观测面相干因子（书图 9.10 复现 + 去倾斜对比）', fontsize=11)
    ax = axs[1]
    r0m = res['r_ax'] / res['r0_mid']
    sel = r0m <= 12
    ax.plot(r0m[sel], 6.88 * r0m[sel] ** (5.0 / 3), 'k-', lw=2,
            label=r'理论 $6.88(r/r_{0,mid})^{5/3}$（式 9.44）')
    ax.plot(r0m[sel], res['Dmid'][sel], 'k-.', lw=1.5,
            label='仿真的中程屏平均结构函数')
    ax.set_xlabel(r'$|\Delta r| / r_{0,mid}$')
    ax.set_ylabel(r'$D_\phi(|\Delta r|)$ [rad$^2$]')
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(r'(b) 中程屏结构函数（$r_{0,mid}$=%.3f m，斜率 %.2f $\approx$ 5/3）'
                 % (res['r0_mid'], res['slope_mid']), fontsize=11)
    fig.suptitle('图 9-14（自绘仿真）湍流仿真验证（%s：50 km, $C_n^2$=1e-16, 20 次实现）'
                 % res['tag'], fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fp = os.path.join(outdir, 'out_coh_%s.png' % res['tag'])
    fig.savefig(fp, dpi=150); plt.close(fig)
    print('    已保存 %s' % fp)


# =====================================================================
# 工程 A：FT 屏 vs 次谐波屏的结构函数（书图 9.3 复现）
# =====================================================================
def run_screens(outdir, nscr_real=40, N=256, r0=0.1, delta=None,
                L0=100.0, l0=0.01, seed=9):
    os.makedirs(outdir, exist_ok=True)
    if delta is None:
        delta = 2.0 / N                              # Listing 9.4：屏边长 2 m
    rng = np.random.default_rng(seed)
    r1d = np.arange(1, N // 2) * delta
    Dth = 6.88 * (r1d / r0) ** (5.0 / 3)
    acc_ft = np.zeros(N // 2 - 1); acc_sh = np.zeros(N // 2 - 1)
    for i in range(nscr_real):
        phz = ft_phase_screen(r0, N, delta, L0, l0, rng)
        Dq, _ = str_fcn2_ft(phz, delta)
        acc_ft += Dq[0, 1:N // 2]
        lo, hi = ft_sh_phase_screen(r0, N, delta, L0, l0, rng)
        Dq, _ = str_fcn2_ft(lo + hi, delta)
        acc_sh += Dq[0, 1:N // 2]
        if (i + 1) % 10 == 0:
            print('    screens %d/%d' % (i + 1, nscr_real))
    acc_ft /= nscr_real; acc_sh /= nscr_real
    # 大间隔偏差与斜率
    r0m = r1d / r0
    sel = (r0m >= 2) & (r0m <= 8)
    cf_sh = np.polyfit(np.log(r0m[sel]), np.log(acc_sh[sel]), 1)
    sel12 = (r0m >= 8) & (r0m <= 12)
    dev_ft = np.max(np.abs(acc_ft[sel12] - Dth[sel12]) / Dth[sel12])
    dev_sh = np.max(np.abs(acc_sh[sel12] - Dth[sel12]) / Dth[sel12])
    i_r0 = int(round(r0 / delta)) - 1
    print('===== 工程 A：FT 屏 vs 次谐波屏（书图 9.3）=====')
    print('    参数：N=%d, delta=%.4f m, r0=%.2f m, L0=%.0f m, l0=%.2f m, %d 块屏'
          % (N, delta, r0, L0, l0, nscr_real))
    print('    D(r=r0)：理论 6.88；SH = %.2f；FT = %.2f'
          % (acc_sh[i_r0], acc_ft[i_r0]))
    print('    2-8 r0 段 log-log 斜率：SH = %.3f（理论 5/3 = 1.667）' % cf_sh[0])
    print('    大间隔（8-12 r0）最大相对偏差：SH = %.1f%%，FT = %.1f%%'
          % (dev_sh * 100, dev_ft * 100))
    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.plot(r0m, Dth, 'k-', lw=2, label=r'理论 $6.88(r/r_0)^{5/3}$（式 9.44）')
    ax.plot(r0m, acc_sh, 'k-.', lw=1.6, label='仿真：FT + 次谐波（Listing 9.3）')
    ax.plot(r0m, acc_ft, 'k:', lw=1.8, label='仿真：纯 FT（Listing 9.2）')
    ax.set_xlabel(r'$|\Delta r| / r_0$')
    ax.set_ylabel(r'$D_\phi(|\Delta r|)$ [rad$^2$]')
    ax.grid(alpha=0.3); ax.legend(fontsize=9)
    ax.set_title('图 9-11（自绘仿真）相位屏结构函数验证（书 9.3 节图 9.3 复现）\n'
                 '纯 FT 屏低频欠采样 -> 大间隔偏低；次谐波法补回倾斜等低阶模式',
                 fontsize=11)
    fig.tight_layout()
    fp = os.path.join(outdir, 'out_fig93_ft_sh.png')
    fig.savefig(fp, dpi=150); plt.close(fig)
    print('    已保存 %s' % fp)
    np.savez_compressed(os.path.join(outdir, 'out_screensA.npz'),
                        r0m=r0m, Dth=Dth, acc_ft=acc_ft, acc_sh=acc_sh,
                        slope_sh=cf_sh[0], dev_ft=dev_ft, dev_sh=dev_sh,
                        D_r0_ft=acc_ft[i_r0], D_r0_sh=acc_sh[i_r0])


# =====================================================================
# 闪烁趋势汇总：三个 Rytov 数下的测量 vs 预测
# =====================================================================
def run_trend(outdir):
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for tag, lab in [('weak', r'50 km, Cn2=1e-17'),
                     ('p7', r'100 km, Cn2=1e-17'),
                     ('main', r'50 km, Cn2=1e-16')]:
        fp = os.path.join(outdir, 'out_%s.npz' % tag)
        if not os.path.exists(fp):
            print('    缺 %s，跳过（先跑对应命令）' % fp)
            continue
        d = np.load(fp)
        rows.append((tag, lab, float(d['rytov']), float(d['var_chi']),
                     float(d['sigI'])))
    print('===== 闪烁 vs Rytov 方差（campaign 汇总）=====')
    print('    %-22s %10s %12s %12s %10s'
          % ('案例', '预测s2_chi', '实测s2_chi', '实测s2_I', '4*预测'))
    for tag, lab, ry, vc, si in rows:
        print('    %-22s %10.4f %12.4f %12.4f %10.4f'
              % (lab, ry, vc, si, 4 * ry))
    if not rows:
        return
    fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.7))
    x = np.arange(len(rows))
    ry = np.array([r[2] for r in rows])
    vc = np.array([r[3] for r in rows])
    si = np.array([r[4] for r in rows])
    axs[0].plot(x, ry, 'k-o', label=r'预测 $\sigma_\chi^2$（式 9.64）')
    axs[0].plot(x, vc, 'r--s', label=r'仿真实测 $\sigma_\chi^2$')
    axs[0].set_xticks(x)
    axs[0].set_xticklabels([r[1] for r in rows], fontsize=9)
    axs[0].set_yscale('log'); axs[0].grid(alpha=0.3, which='both')
    axs[0].yaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, p: '%g' % v))    # 避免 mathtext 负号缺字形
    axs[0].set_ylabel(r'$\sigma_\chi^2$（对数振幅方差）')
    axs[0].legend(fontsize=9)
    axs[0].set_title('(a) 对数振幅方差随 Rytov 数', fontsize=11)
    axs[1].plot(x, 4 * ry, 'k-o', label=r'$4\sigma_\chi^2$（弱起伏 $\sigma_I^2$）')
    axs[1].plot(x, si, 'r--s', label=r'仿真实测 $\sigma_I^2$（闪烁指数）')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels([r[1] for r in rows], fontsize=9)
    axs[1].set_yscale('log'); axs[1].grid(alpha=0.3, which='both')
    axs[1].yaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, p: '%g' % v))
    axs[1].set_ylabel(r'$\sigma_I^2$（闪烁指数）'); axs[1].legend(fontsize=9)
    axs[1].set_title('(b) 闪烁指数：弱起伏区吻合，强起伏区饱和', fontsize=11)
    fig.suptitle('图 9-15（自绘仿真）闪烁统计随 Rytov 方差的变化'
                 '（每例 20 次实现；口径内平均）', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fp = os.path.join(outdir, 'out_scint_trend.png')
    fig.savefig(fp, dpi=150); plt.close(fig)
    print('    已保存 %s' % fp)


def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'main'
    outdir = argv[2] if len(argv) > 2 else _HERE
    if cmd == 'screens':
        run_screens(outdir)
    elif cmd == 'main':
        run_campaign('main', outdir, D2=0.5, Dz=50e3, Cn2=1e-16,
                     nscr=11, nreals=20, seed=9)
    elif cmd == 'p7':
        run_campaign('p7', outdir, D2=0.5, Dz=100e3, Cn2=1e-17,
                     nscr=7, nreals=20, seed=11)
    elif cmd == 'weak':
        run_campaign('weak', outdir, D2=0.5, Dz=50e3, Cn2=1e-17,
                     nscr=7, nreals=20, seed=12)
    elif cmd == 'trend':
        run_trend(outdir)
    else:
        print(__doc__)


if __name__ == '__main__':
    main(sys.argv)
