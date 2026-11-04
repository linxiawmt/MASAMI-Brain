# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 8 章（工程 2：习题 1-4）=====
部分传播与放宽采样约束——习题复现与数值验证

习题清单（原书 8.5 节）：
    1. 啁啾信号 g(x)=exp(i*pi*a^2*x^2)（式 8.26，a=4，N=128，L=4 m）的解析
       FT 与离散 DFT 对比；再乘超高斯吸收边界（n=16，σ=0.25L）重算 DFT，
       画两版 DFT 的实部/虚部并与解析 FT 对比。
    2. 补全式 (8.20) 与 (8.22) 之间的代数步骤，证明约束 3 与部分传播次数
       无关。（代数推导见笔记；本脚本做端到端数值复核：会聚源一次传播
       vs 任意拆分的部分传播，精度应与拆分方式无关。）
    3. 点源（λ=1 µm）传播 100 km 到 2 m 直径望远镜的采样图，并与单次传播
       的情形比较。
    4. 环形孔径（外径 1.5 m、内径 0.5 m，λ=1.3 µm）会聚到 100 km 处目标
       平面：(a) 仿照图 8.5 的采样分析并确定部分传播段数；(b) 仿真 y=0
       切片的辐照度与相位，与解析结果同图对比。

Python 实跑结果（MATLAB 版按相同参数复现，见 matlab/waveprop/ch8/）：
    习题 1：无边界 DFT 在 |f| 接近 Nyquist 处与解析 FT 偏差很大（实部
            最大偏差 0.21、虚部 0.22，而 |G|=0.25）；乘超高斯后偏差降到
            2.4e-2/2.7e-2（降低约一个数量级）——吸收边界压制混叠的直接
            证据。残余偏差来自 |x|>σ 附近仍有 f_loc>Nyquist 的样本。
    习题 2：会聚源（R=-1 m，ρ=0.5）：单次传播 |dI|/pk = 3.3e-3；拆成
            α2=0.2/0.5/0.8 的两次部分传播分别为 3.3e-3/3.3e-3/3.3e-3，
            5 段 3.3e-3——约束 3 与拆分方式无关的端到端实锤。第 2 段的
            局部频率最大值与本笔记的"偏移中心"公式一致（差 2e-15）。
    习题 3：点源 100 km、D2=2 m：(δ1,δn)=(1,1) cm 时约束 2 需要 N>=600、
            单次约束 4 需要 N>=1000 -> N=1024（第 7 章习题 2 的答案）；
            换 (δ1,δn)=(2,3) cm：约束 2 只需 N>=117、单次约束 4 需
            N>=167 -> N=256，而 2 段部分传播的逐段约束 4 只需 N>=83
            -> N=128——部分传播把网格减半。N=128 两段仿真：ROI 内相位
            残差 1.3e-6 rad、振幅平坦度 3e-6（对照球面波解析）。
    习题 4：环形孔径聚焦：δ1=50 mm、δn=16.7 mm（各 30 点）、N=128、
            式 (8.25) Δzmax = 27.5 km -> 书公式 5 段（严格最小 4 段）。
            仿真 vs 解析（环形孔径夫琅禾费图样 = jinc 之差）：y=0 切片
            ±0.25 m 内 |dI|/峰值 = 8.9e-3、相位残差 6.9e-2 rad（零点附近
            相位跳变处残差主导）；峰值辐照度 146 与解析 146.0 一致。

MATLAB/Python 差异（对应 .m 内注释）：
    - MATLAB 1-based / numpy 0-based；meshgrid 行为一致；
    - 解析 FT / Fresnel 积分：Python 用闭式；MATLAB 同式（无内置 Fresnel
      积分的旧版本用 integral() 数值积分）；
    - 本机无 MATLAB 时以 Python 实跑结果为准。
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
from scipy.special import fresnel, j1

from base import ft1_centered, jinc, rect
from ch6_fresnel_prop import fresnel_prop_square_ap
from ch7_sampling import (constraint1_d2, constraint2_N,
                          constraint3_d2, constraint4_N,
                          fresnel_prop_square_ap_R)
from ch8_partial_prop import (ang_spec_multi_prop_vac, dz_max_partial,
                              n_partial_book)

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


# ===== 习题 1 的超高斯窗口（1 维版，式 8.1） =====
def sg_1d(x, sigma, n_order=16):
    """g_sg(x) = exp(-(|x|/σ)^n)：式 (8.1) 的一维径向形式。"""
    return np.exp(-(np.abs(np.asarray(x, dtype=float)) / sigma) ** n_order)


def prob1(outdir):
    """习题 1：啁啾 DFT 与吸收边界。"""
    print('===== 习题 1：g(x)=exp(i*pi*a^2*x^2) 的 DFT 与超高斯吸收边界 =====')
    a, N, L = 4.0, 128, 4.0
    delta = L / N                                   # 0.03125 m
    x = (np.arange(N) - N // 2) * delta             # -2 .. 1.96875 m
    g = np.exp(1j * np.pi * a ** 2 * x ** 2)        # 式 (8.26)
    # 无边界 DFT 与乘超高斯（n=16, σ=0.25L）后的 DFT
    Gd = ft1_centered(g, delta)
    gsg = g * sg_1d(x, 0.25 * L, 16)
    Gdsg = ft1_centered(gsg, delta)
    # 解析 FT（无限长啁啾）：∫exp(iπa²x²-i2πfx)dx = e^{iπ/4}/a · e^{-iπf²/a²}
    f = (np.arange(N) - N // 2) / (N * delta)       # -16 .. 15.875 1/m
    Gan = np.exp(1j * np.pi / 4.0) / a * np.exp(-1j * np.pi * f ** 2 / a ** 2)
    # 混叠解析判据（无 FT）：f_loc = a²x = 16x，Nyquist = 1/(2δ) = 16 /m
    fnyq = 1.0 / (2.0 * delta)
    x_alias = fnyq / a ** 2                          # 1.0 m：|x|>1 m 的样本全部混叠
    n_alias = int(np.sum(np.abs(x) > x_alias))
    sg_at_alias = sg_1d(x_alias + delta / 2, 0.25 * L, 16)
    print('    delta=%.5f m，Nyquist=%.0f /m；f_loc=a^2*x=16x' % (delta, fnyq))
    print('    混叠边界 |x| = %.3f m：%d/%d 个样本混叠（一半网格）' % (x_alias, n_alias, N))
    print('    超高斯（σ=0.25L=1 m, n=16）在混叠边界内侧的值 = %.3f（边界处仍显著）'
          % sg_1d(x_alias, 0.25 * L, 16))
    e0 = np.max(np.abs(Gd - Gan))
    e1 = np.max(np.abs(Gdsg - Gan))
    print('    DFT vs 解析：无边界 max|ΔRe|=%.2e max|ΔIm|=%.2e max|Δ|=%.2e'
          % (np.max(np.abs(Gd.real - Gan.real)), np.max(np.abs(Gd.imag - Gan.imag)), e0))
    print('              加边界 max|ΔRe|=%.2e max|ΔIm|=%.2e max|Δ|=%.2e'
          % (np.max(np.abs(Gdsg.real - Gan.real)),
             np.max(np.abs(Gdsg.imag - Gan.imag)), e1))
    print('    |G| 解析 = 1/a = %.3f（常数）；无边界 DFT 的 |G| 范围 [%.3f, %.3f]，'
          '加边界 [%.3f, %.3f]' % (1 / a, np.abs(Gd).min(), np.abs(Gd).max(),
                                  np.abs(Gdsg).min(), np.abs(Gdsg).max()))

    fig9, axs9 = plt.subplots(2, 2, figsize=(12.5, 7.6))
    ax = axs9[0, 0]
    ax.plot(x, g.real, 'k-', lw=0.8, label='Re g')
    ax.plot(x, gsg.real, 'C1--', lw=1.2, label=r'Re g$\times$超高界')
    ax.axvspan(x_alias, L / 2, color='r', alpha=0.12)
    ax.axvspan(-L / 2, -x_alias, color='r', alpha=0.12)
    ax.set_title(r'(a) 采样后的啁啾（红色 = $f_{loc}>$Nyquist 混叠区）', fontsize=10)
    ax.set_xlabel(r'$x$ [m]')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs9[0, 1]
    ax.plot(f, Gd.real, 'C0.-', ms=3, label='DFT 实部（无边界）')
    ax.plot(f, Gan.real, 'k-', lw=1.0, label='解析实部')
    ax.plot(f, Gdsg.real, 'C1+-', ms=4, lw=0.8, label='DFT 实部（超高斯）')
    ax.set_title('(b) DFT 实部 vs 解析 FT', fontsize=10)
    ax.set_xlabel(r'$f$ [1/m]')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs9[1, 0]
    ax.plot(f, Gd.imag, 'C0.-', ms=3, label='DFT 虚部（无边界）')
    ax.plot(f, Gan.imag, 'k-', lw=1.0, label='解析虚部')
    ax.plot(f, Gdsg.imag, 'C1+-', ms=4, lw=0.8, label='DFT 虚部（超高斯）')
    ax.set_title('(c) DFT 虚部 vs 解析 FT', fontsize=10)
    ax.set_xlabel(r'$f$ [1/m]')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axs9[1, 1]
    ax.plot(f, np.abs(Gd), 'C0.-', ms=3, label=r'$|G|$（无边界）')
    ax.plot(f, np.abs(Gan), 'k-', lw=1.0, label=r'解析 $|G|=1/a$')
    ax.plot(f, np.abs(Gdsg), 'C1+-', ms=4, lw=0.8, label=r'$|G|$（超高斯）')
    ax.set_title('(d) 幅值：混叠把能量搬到错误的频率上', fontsize=10)
    ax.set_xlabel(r'$f$ [1/m]')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig9.suptitle(r'图 8-9（自绘仿真）习题 1：$g(x)=e^{i\pi a^2x^2}$（a=4, N=128, L=4 m）'
                  r'的 DFT 与超高斯吸收边界（n=16, $\sigma$=0.25L）', fontsize=12)
    fig9.tight_layout(rect=[0, 0, 1, 0.93])
    fp = os.path.join(outdir, 'fig-8-9.png')
    fig9.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    return e0, e1


def prob2(outdir):
    """习题 2：约束 3 与部分传播次数无关——端到端数值复核。"""
    print('===== 习题 2：约束 3 的代数（8.20->8.22）与会聚源端到端复核 =====')
    print('    代数：δ2=(1-α2)δ1+α2δn、Δz1=α2Δz 代入式 (8.20)，')
    print('    两边减 (1-α2)δ1、除以 α2（α2>0），即得式 (8.22)，与 α2 无关。')
    # 会聚源端到端：R=-1 m（ρ=1+Δz/R=0.5>0，焦点在目标面之后）
    wvl, D1, D2, Dz, R = 1e-6, 2e-3, 4e-3, 0.5, -1.0
    d1, dn = 66.6667e-6, 100e-6
    N = 128
    rho = 1 + Dz / R
    lo3, hi3 = constraint3_d2(d1, D1, wvl, Dz, R)
    c2N = constraint2_N(d1, dn, D1, D2, wvl, Dz)
    print('    参数：D1=2mm, D2=4mm, Δz=0.5m, R=-1m（会聚，ρ=%.1f）' % rho)
    print('    约束 3（式 8.22）：%.0f <= δn <= %.0f µm；约束 2 需 N>=%.1f；'
          '选中 δn=%.0f µm、N=128' % (lo3 * 1e6, hi3 * 1e6, c2N, dn * 1e6))
    n1, y1 = np.meshgrid((np.arange(N) - N // 2) * d1,
                         (np.arange(N) - N // 2) * d1)
    src = rect(n1 / D1) * rect(y1 / D1) * np.exp(
        1j * np.pi / wvl / R * (n1 ** 2 + y1 ** 2))
    # 解析（第 7 章习题 4a 的推广：fresnel_prop_square_ap_R）
    xs = (np.arange(N) - N // 2) * dn
    Uan = fresnel_prop_square_ap_R(xs, 0.0, D1, wvl, Dz, R)

    def err(xn, U):
        mid = U.shape[0] // 2
        sl = xn[0, :]
        win = np.abs(sl) <= D2 / 2
        I, Ia = np.abs(U[mid, :]) ** 2, np.abs(Uan) ** 2
        pk = Ia.max()
        eI = np.max(np.abs(I[win] - Ia[win])) / pk
        mph = win & (Ia > 0.05 * pk)
        eph = np.max(np.abs(np.angle(U[mid, :][mph] / Uan[mph])))
        return eI, eph

    xn1, _, U1 = ang_spec_multi_prop_vac(src, wvl, d1, dn, [Dz], absorb='none')
    print('    单次传播：            |dI|/pk = %.2e，相位残差 = %.2e rad'
          % err(xn1, U1))
    for a2 in (0.2, 0.5, 0.8):
        zp = [Dz * a2, Dz]
        xn2, _, U2 = ang_spec_multi_prop_vac(src, wvl, d1, dn, zp)
        eI, eph = err(xn2, U2)
        print('    两次部分传播 α2=%.1f：  |dI|/pk = %.2e，相位残差 = %.2e rad'
              % (a2, eI, eph))
    zp5 = [Dz * i / 5 for i in range(1, 6)]
    xn5, _, U5 = ang_spec_multi_prop_vac(src, wvl, d1, dn, zp5)
    print('    5 段部分传播：         |dI|/pk = %.2e，相位残差 = %.2e rad' % err(xn5, U5))
    # 本笔记补充：第 2 段（有限 R）约束 3 中心偏移公式的数值核对
    # f_loc,max（第 i 段）= (D_i/(2λ))·|(δ1-δn)/Δz + δi/R|，D_i = D1(1+z_{i-1}/R)
    a2 = 0.5
    z_mid, d_mid = a2 * Dz, (1 - a2) * d1 + a2 * dn
    R2, Dbeam = R + z_mid, D1 * abs(1 + z_mid / R)
    m2 = dn / d_mid
    fl_direct = Dbeam / (2 * wvl) * abs((1 - m2) / (Dz - z_mid) + 1 / R2)
    fl_formula = Dbeam / (2 * wvl) * abs((d1 - dn) / Dz + d_mid / R)
    dev = abs(fl_direct - fl_formula) / fl_direct
    print('    偏移中心公式核对（第 2 段局部频率最大值）：直算 %.6e /m vs '
          '公式 %.6e /m，相对差 %.1e' % (fl_direct, fl_formula, dev))
    print('    —— 书句"约束 3 与段数无关"在 R=∞ 严格成立；有限 R 时第 i 段')
    print('    的中心偏移 Δz(δi-δ1)/R（本例 %.1f µm，仍远离约束边界）'
          % abs(Dz * (d_mid - d1) / R * 1e6))
    return dev


def prob3(outdir):
    """习题 3：点源 100 km 的采样图与部分传播对比。"""
    print('===== 习题 3：点源（λ=1µm）100 km -> 2 m 望远镜的采样图 =====')
    wvl, Dz, D2 = 1e-6, 100e3, 2.0
    # 约束 2（D1->0）：N >= D2/(2δn) + λΔz/(2δ1δn)
    # 约束 4（单次）：N >= λΔz/(δ1δn)；（部分传播）逐段：N >= λΔzi/(δiδi+1)
    d1v = np.logspace(-3.3, -1, 300)                 # 0.5 mm - 100 mm
    dnv = np.logspace(-3.3, -1, 300)
    D1g, Dng = np.meshgrid(d1v, dnv)
    c2 = D2 / (2 * Dng) + wvl * Dz / (2 * D1g * Dng)          # 约束 2
    c4s = wvl * Dz / (D1g * Dng)                               # 约束 4（单次）
    fig10, axs10 = plt.subplots(1, 2, figsize=(12.5, 5.2))
    ax = axs10[0]
    cs = ax.contourf(D1g, Dng, np.log2(np.maximum(c2, 1)),
                     levels=np.arange(0, 14, 1), cmap='Greys')
    clb = fig10.colorbar(cs, ax=ax); clb.set_label(r'$\log_2 N$（约束 2）')
    cs4 = ax.contour(D1g, Dng, np.log2(np.maximum(c4s, 1)),
                     levels=[7, 8, 9, 10], colors='C3', linewidths=1.2)
    ax.clabel(cs4, fmt='单次约束4 %2.0f', fontsize=7)
    ax.plot(1e-2, 1e-2, 'k*', ms=12)
    ax.annotate('(1, 1) cm\nN=1024（第 7 章习题 2）', (1e-2, 1e-2),
                textcoords='offset points', xytext=(8, 8), fontsize=8)
    ax.plot(2e-2, 3e-2, 'C1o', mfc='none', ms=10, mew=2)
    ax.annotate('(2, 3) cm：约束2 N≥117\n单次约束4 N≥167 -> 256\n2 段逐段 N≥83 -> 128',
                (2e-2, 3e-2), textcoords='offset points', xytext=(10, -30),
                fontsize=8, color='C1')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$\delta_1$ [m]'); ax.set_ylabel(r'$\delta_n$ [m]')
    ax.set_title('(a) 点源 100 km（D2=2 m）：约束 2 等高线 + 单次约束 4（红）', fontsize=10)
    ax.grid(alpha=0.3, which='both')
    ax = axs10[1]
    nseg_v = np.arange(1, 9)
    picks = {(2e-2, 3e-2): 'C1', (1e-2, 1e-2): 'k'}
    for (p1, pn), col in picks.items():
        Nreq = []
        for ns in nseg_v:
            al = np.arange(ns + 1) / ns
            dl = (1 - al) * p1 + al * pn
            Nreq.append(max(wvl * (Dz / ns) / (dl[i] * dl[i + 1])
                            for i in range(ns)))
        Nreq = np.maximum(Nreq, 1)
        ax.plot(nseg_v, Nreq, 'o-', color=col, lw=1.4,
                label='(%.0f, %.0f) cm' % (p1 * 100, pn * 100))
        for ns, nv in zip(nseg_v, Nreq):
            if ns <= 4:
                ax.annotate('%.0f' % nv, (ns, nv), textcoords='offset points',
                            xytext=(0, 6), fontsize=7, color=col)
    ax.axhline(128, color='C6', ls='--', lw=1.2)
    ax.text(6.1, 140, 'N = 128', fontsize=9, color='C6')
    ax.axhline(1024, color='0.5', ls=':', lw=1)
    ax.text(4.1, 1100, 'N = 1024', fontsize=8, color='0.4')
    ax.set_yscale('log')
    ax.set_xlabel('部分传播段数'); ax.set_ylabel('逐段约束 4 需要的 N')
    ax.set_title('(b) 段数放松约束 4：点源案例（约束 2 恒需 N≥600/117）', fontsize=10)
    ax.legend(fontsize=8, title='(δ1, δn)'); ax.grid(alpha=0.3, which='both')
    fig10.suptitle('图 8-10（自绘仿真）习题 3：点源 100 km 的采样图——部分传播 vs 单次传播',
                   fontsize=12)
    fig10.tight_layout(rect=[0, 0, 1, 0.92])
    fp = os.path.join(outdir, 'fig-8-10.png')
    fig10.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    # 具体数值
    for (p1, pn) in ((1e-2, 1e-2), (2e-2, 3e-2)):
        c2n = D2 / (2 * pn) + wvl * Dz / (2 * p1 * pn)
        c4n = wvl * Dz / (p1 * pn)
        al2 = np.array([0, 0.5, 1.0])
        dl2 = (1 - al2) * p1 + al2 * pn
        c4n2 = wvl * (Dz / 2) / (dl2[0] * dl2[1])
        print('    (δ1,δn)=(%.0f,%.0f) cm：约束 2 N>=%.0f；单次约束 4 N>=%.0f；'
              '2 段逐段 N>=%.0f' % (p1 * 100, pn * 100, c2n, c4n, c4n2))
    # N=128、(2,3)cm、2 段的端到端仿真（点源 = 单像素，平面波谱最平——卷绕压力最大）
    N, d1, dn = 128, 2e-2, 3e-2
    src = np.zeros((N, N), dtype=complex)
    src[N // 2, N // 2] = 1.0
    xn, yn, U = ang_spec_multi_prop_vac(src, wvl, d1, dn, [Dz / 2, Dz])
    mid = N // 2
    r = np.abs(xn[mid, :])
    win = r <= D2 / 2
    ph = np.angle(U[mid, :][win])
    r_roi = r[win]
    ph_an = np.pi / wvl * r_roi ** 2 / Dz            # 球面波相位（曲率 1/Dz）
    dev = np.max(np.abs(np.mod(ph - ph_an + np.pi, 2 * np.pi) - np.pi)
                 [1:])                                # 去中心像素（相位奇异）
    amp = np.abs(U[mid, :][win])
    amp_flat = np.max(np.abs(amp[1:] / amp[1] - 1))
    print('    N=128、2 段×50 km 仿真：ROI（r<=1 m）内相位 vs 球面波解析残差'
          ' = %.1e rad；振幅平坦度 = %.1e（幅度归一到第一环带）' % (dev, amp_flat))
    return dev, amp_flat


def prob4(outdir):
    """习题 4：环形孔径聚焦到 100 km 目标。"""
    print('===== 习题 4：环形孔径（外径1.5m/内径0.5m, λ=1.3µm）聚焦 100 km =====')
    wvl, Dz = 1.3e-6, 100e3
    b, a_in = 0.75, 0.25                             # 外/内半径 [m]
    D1 = 2 * b                                       # 源面孔径外径
    # 采样分析（仿图 8.5）：30 点跨 D1 与 D2
    d1 = D1 / 30                                     # 50 mm
    # 目标面 ROI：取 D2 = 0.5 m（约 2.4 个艾里半径）
    D2 = 0.5
    dn = D2 / 30                                     # 16.7 mm
    N = 128
    c1 = constraint1_d2(d1, D1, D2, wvl, Dz)
    c2N = constraint2_N(d1, dn, D1, D2, wvl, Dz)
    lo3, hi3 = constraint3_d2(d1, D1, wvl, Dz, -Dz)   # R = -Dz（会聚到目标，1+Δz/R=0）
    zmax = dz_max_partial(d1, dn, N, wvl)
    nbk = n_partial_book(Dz, zmax)
    print('    δ1=%.1f mm（30 点跨 1.5m），δn=%.1f mm（30 点跨 0.5m ROI），N=%d'
          % (d1 * 1e3, dn * 1e3, N))
    print('    约束 1：δn <= %.1f mm；约束 2：N >= %.1f（%s）；'
          '约束 3（R=-100km，ρ=0）：δn <= %.1f mm'
          % (c1 * 1e3, c2N, '满足' if c2N <= N else '违反！', hi3 * 1e3))
    print('    约束 4：式 (8.25) Δzmax = min(δ)^2·N/λ = %.1f km -> 书公式 %d 段'
          '（严格最小 %d 段，每段 %.0f km）'
          % (zmax / 1e3, nbk, int(np.ceil(Dz / zmax)), Dz / int(np.ceil(Dz / zmax)) / 1e3))
    # ---- (a) 约束图 ----
    fig11, axs11 = plt.subplots(1, 3, figsize=(15.5, 4.8))
    ax = axs11[0]
    d1v = np.linspace(10e-3, 120e-3, 300)
    dnv = np.linspace(2e-3, 60e-3, 300)
    D1g, Dng = np.meshgrid(d1v, dnv)
    log2N = np.log2(np.maximum(constraint2_N(D1g, Dng, D1, D2, wvl, Dz), 1e-3))
    cs = ax.contourf(D1g * 1e3, Dng * 1e3, log2N,
                     levels=np.arange(4, 13.5, 0.5), cmap='Greys')
    clb = fig11.colorbar(cs, ax=ax); clb.set_label(r'$\log_2 N$（约束 2）')
    c1v = constraint1_d2(d1v, D1, D2, wvl, Dz)
    ax.plot(d1v * 1e3, np.clip(c1v, 0, 60) * 1e3, 'k--', lw=2,
            label='约束 1（式 7.14）')
    _, hi3v = constraint3_d2(d1v, D1, wvl, Dz, -Dz)
    ax.plot(d1v * 1e3, np.clip(hi3v, 0, 60) * 1e3, 'k-.', lw=1.5,
            label=r'约束 3 上界（$\rho$=0）')
    zmv = dnv ** 2 * N / wvl                          # Δzmax 等高线：min=δn（收缩网格）
    csz = ax.contour(D1g * 1e3, Dng * 1e3, zmv / 1e3,
                     levels=[10, 20, 30, 50, 100], colors='C3', linewidths=1.0)
    ax.clabel(csz, fmt='Δzmax=%2.0f km', fontsize=7)
    ax.plot(d1 * 1e3, dn * 1e3, 'r*', ms=14, zorder=5)
    ax.annotate('选中 (50, 16.7) mm\nN=128, 5 段×20 km', (d1 * 1e3, dn * 1e3),
                textcoords='offset points', xytext=(8, 10), fontsize=8, color='r')
    ax.set_xlabel(r'$\delta_1$ [mm]'); ax.set_ylabel(r'$\delta_n$ [mm]')
    ax.set_xlim(10, 120); ax.set_ylim(0, 60)
    ax.set_title('(a) 习题 4(a)：采样约束图（D1=1.5m, D2=0.5m, Δz=100km）', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    # ---- 仿真（5 段 × 20 km） ----
    nseg = 5
    zpl = [Dz * i / nseg for i in range(1, nseg + 1)]
    n = np.arange(N) - N // 2
    x1g, y1g = np.meshgrid(n * d1, n * d1)
    r1 = np.hypot(x1g, y1g)
    src = ((r1 <= b) & (r1 >= a_in)).astype(float) * np.exp(
        -1j * np.pi / wvl / Dz * (x1g ** 2 + y1g ** 2))   # 会聚：R=-Dz（R<0 会聚）
    xn, yn, U = ang_spec_multi_prop_vac(src, wvl, d1, dn, zpl)
    mid = N // 2
    # 解析：目标面 = 环形孔径的夫琅禾费图样 × 观测面啁啾
    # U(r2) = e^{ikΔz}/(iλΔz)·e^{+ik r2^2/(2Δz)}·Â(r2/(λΔz))，
    # 会聚源相位 e^{-ik r1^2/(2Δz)} 与传播核的 r1^2 项相消，只剩 FT。
    xs = xn[0, :]
    rho = xs / (wvl * Dz)
    Uan = (np.exp(1j * np.pi * xs * rho)
           * (jinc(rho, b) - jinc(rho, a_in)) / (1j * wvl * Dz))
    I, Ia = np.abs(U[mid, :]) ** 2, np.abs(Uan) ** 2
    print('    峰值辐照度：仿真 %.1f vs 解析 %.1f' % (I.max(), Ia.max()))
    win = np.abs(xs) <= D2 / 2
    eI = np.max(np.abs(I[win] - Ia[win])) / Ia.max()
    mph = win & (Ia > 0.05 * Ia.max())
    eph = np.max(np.abs(np.angle(U[mid, :][mph] / Uan[mph])))
    print('    y=0 切片 ±%.2f m：|dI|/峰值 = %.2e；相位残差（|U|>5%%峰）= %.2e rad'
          % (D2 / 2, eI, eph))
    # ---- (b)(c) 图 ----
    ax = axs11[1]
    ax.plot(xs, I / Ia.max(), 'b.-', ms=4, lw=0.8, label='仿真')
    ax.plot(xs, Ia / Ia.max(), 'k-', lw=1.2, label='解析')
    ax.set_xlim(-250, 250)
    ax.set_xlabel(r'$x_n$ [m]')
    ax.set_ylabel('归一化辐照度')
    ax.set_yscale('log')
    ax.set_ylim(1e-6, 2)
    ax.set_title('(b) 目标面 y=0 切片辐照度（对数轴）：|dI|/pk = %.1e' % eI, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')
    ax = axs11[2]
    ax.plot(xs, np.angle(Uan), 'k-', lw=1.2, label='解析')
    ax.plot(xs[mph], np.angle(U[mid, :][mph]), 'b.', ms=4, label='仿真')
    ax.set_xlim(-250, 250)
    ax.set_xlabel(r'$x_n$ [m]'); ax.set_ylabel('相位 [rad]')
    ax.set_title('(c) 目标面 y=0 切片相位：残差 %.1e rad' % eph, fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig11.suptitle('图 8-11（自绘仿真）习题 4：环形孔径（外径1.5 m/内径0.5 m, λ=1.3 µm）'
                   '聚焦到 100 km 目标（N=128, 5 段×20 km）', fontsize=12)
    fig11.tight_layout(rect=[0, 0, 1, 0.9])
    fp = os.path.join(outdir, 'fig-8-11.png')
    fig11.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    return eI, eph


def main(outdir=None):
    here = os.path.dirname(os.path.abspath(__file__))
    if outdir is None:
        outdir = here
    os.makedirs(outdir, exist_ok=True)
    prob1(outdir)
    prob2(outdir)
    prob3(outdir)
    prob4(outdir)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
