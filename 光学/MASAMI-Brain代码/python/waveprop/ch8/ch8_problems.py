# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 8 章（工程 2：习题 1-4）=====
部分传播与放宽采样约束——习题复现与数值验证（原书 8.5 节）。

习题清单：
    1. 啁啾信号 g(x)=exp(i*pi*a^2*x^2)（式 8.26，a=4，N=128，L=4 m）的解析
       FT 与离散 DFT 对比；再把信号预乘 super-Gaussian 吸收边界
       （n=16、sigma=0.25L）重算 DFT，画两版 DFT 的实部/虚部并与解析 FT 对比。
       要点：a=4、L=4 m 时 f_loc = a^2*x = 16x，Nyquist = 16/m，混叠边界
       |x| = 1 m = L/4；sigma = 0.25L = 1 m 的窗恰好在混叠区把信号压掉
       ——吸收边界抑制混叠的一维玩具版。
    2. 补全式 (8.20) -> (8.22) 之间的代数步骤（alpha_2 消去），证明约束 3
       与分段数无关。代数推导见笔记；工程 1 的验证 0 已做 1000 组随机参数
       数值复核，这里再打印推导要点。
    3. 点源（lambda=1 um）传播 100 km 到 2 m 直径望远镜孔径的采样图
       （约束 1&2 等高线 + 可行选点），并与"只有一次传播"的情形比较：
       delta1=delta2=1 cm、N=1024 时 Delta_z_max = min^2*N/lambda = 1.02e5 m
       >= 1e5 m —— 单次传播约束 4 本来就不紧，分段在这里没有收益；
       分段的收益出现在网格间距小（delta ~ 几十 um）的场景（书例）。
    4. 环形孔径（外径 1.5 m、内径 0.5 m、lambda=1.3 um）平面波聚焦到
       100 km 处目标平面：(a) 采样分析（约束 1/2/3 + 分段数确定）；
       (b) yn=0 切片辐照度与相位，仿真 vs 解析同图。
       解析：源带聚焦相位 exp(-i*k*r1^2/(2*Dz)) 时菲涅耳积分的二次相位
       恰好相消，观察面场 = e^{ikDz} e^{ikr^2/(2Dz)} /(i*lambda*Dz) x
       环孔光瞳的傅里叶变换（精确成立，非近似）：
         P(r) = circ(2r/D_out) - circ(2r/D_in)
         P~(f) = (pi/4)[D_out^2 * 2J1(pi D_out f)/(pi D_out f)
                        - D_in^2 * 2J1(pi D_in f)/(pi D_in f)]，f = r/(lambda*Dz)。

数值验证（运行本脚本全部自动打印，图存 out_*.png）：
    (1) 习题 1：无窗 DFT vs 加窗 DFT vs 解析 FT（实部/虚部），量化误差 ->
        out_fig10_prob1.png（笔记 fig-8-10）。
    (2) 习题 2：打印推导 + 引用工程 1 验证 0 的随机参数结论。
    (3) 习题 3：点源采样图 + 单次/分段约束 4 对照 -> out_fig11_prob3.png
        （笔记 fig-8-11）。
    (4) 习题 4：采样分析（约束 2 等高线 + 选点）+ yn=0 辐照度/相位
        仿真 vs 解析 -> out_fig12_prob4.png（笔记 fig-8-12）。

运行：C:/Anocanda/envs/machinelearn/python.exe -I ch8_problems.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import j1

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

from matplotlib.ticker import FuncFormatter


def plain_log(axis):
    """log 轴用纯文本刻度（避免 mathtext 的 U+2212 负号被 SimHei 顶成缺字框）。"""
    axis.set_major_formatter(FuncFormatter(lambda v, p: '%g' % v))
    axis.set_minor_formatter(FuncFormatter(lambda v, p: ''))


_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                # -I 模式下脚本目录不在 path
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'fourier', 'common'))
from base import rect, circle, ft1_centered as ft1       # noqa: E402
from ch8_partial_prop import (                           # noqa: E402
    ang_spec_multi_prop_vac, super_gaussian_absorber,
    constraint1_d2, constraint2_N, constraint3_d2,
    dz_max_partial, delta_middle)


def prob1(outdir):
    """习题 1：啁啾信号的 DFT 与吸收边界预乘。"""
    print('===== 习题 1：g(x)=exp(i*pi*a^2*x^2) 的 DFT 与吸收边界预乘 =====')
    a = 4.0
    N = 128
    L = 4.0
    delta = L / N
    n = np.arange(N) - N // 2
    x = n * delta
    g = np.exp(1j * np.pi * a ** 2 * x ** 2)
    # super-Gaussian 吸收边界（n=16、sigma=0.25L；r 取索引半径 |x|/delta）
    sig_idx = 0.25 * L / delta
    w = np.exp(-(np.abs(x) / delta) ** 16 / sig_idx ** 16)   # 一维超高斯
    g_win = w * g
    # 解析 FT：FT{exp(i pi a^2 x^2)} = e^{i pi/4}/a * exp(-i pi f^2/a^2)
    df = 1.0 / (N * delta)                               # 频域采样间隔（ft1_centered 约定）
    f = n * df
    G_an = np.exp(1j * np.pi / 4.0) / a * np.exp(-1j * np.pi * f ** 2 / a ** 2)
    G_num = ft1(g, delta)
    G_win = ft1(g_win, delta)
    # 量化（比较低频段 |f| <= 8 m^-1，即窗明显起作用的范围）
    band = np.abs(f) <= 8.0
    e_raw = np.max(np.abs(G_num[band] - G_an[band]))
    e_win = np.max(np.abs(G_win[band] - G_an[band]))
    print('    delta=%.4f m, Nyquist=%.1f /m；f_loc = a^2*x = 16x，混叠边界 |x| = %.2f m'
          % (delta, 1 / (2 * delta), 1.0 / (2 * a ** 2 * delta)))
    print('    窗在 |x|=1 m（=混叠边界）处 |w| = %.3f；|x|=2 m 处 %.2e'
          % (np.interp(1.0, x[x > 0], w[x > 0]), np.interp(2.0, x[x > 0], w[x > 0])))
    print('    低频段 |f|<=8 与解析 FT 的最大偏差：无窗 DFT %.3f，加窗 DFT %.3f'
          '（振幅尺度 1/a = %.2f）' % (e_raw, e_win, 1 / a))
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2))
    axs[0].plot(x, w, 'k-', lw=1.2)
    axs[0].set_title('(a) super-Gaussian 窗（n=16, σ=0.25L）', fontsize=10)
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel('g_sg(x)')
    axs[0].axvspan(-1, 1, color='0.9')
    axs[0].text(0, 0.5, '混叠区 |x|>1 m\n被窗压掉', ha='center', fontsize=8)
    axs[0].grid(alpha=0.3)
    axs[1].plot(f, G_an.real, 'k-', lw=1.3, label='解析 Re')
    axs[1].plot(f, G_num.real, 'r.', ms=4, label='无窗 DFT Re')
    axs[1].plot(f, G_win.real, 'b.', ms=3, label='加窗 DFT Re')
    axs[1].set_xlim(-16, 16)
    axs[1].set_xlabel('f [m$^{-1}$]')
    axs[1].set_title('(b) 实部：无窗 DFT 混叠抬起"本应平坦"的谱', fontsize=10)
    axs[1].legend(fontsize=8)
    axs[1].grid(alpha=0.3)
    axs[2].plot(f, G_an.imag, 'k-', lw=1.3, label='解析 Im')
    axs[2].plot(f, G_num.imag, 'r.', ms=4, label='无窗 DFT Im')
    axs[2].plot(f, G_win.imag, 'b.', ms=3, label='加窗 DFT Im')
    axs[2].set_xlim(-16, 16)
    axs[2].set_xlabel('f [m$^{-1}$]')
    axs[2].set_title('(c) 虚部', fontsize=10)
    axs[2].legend(fontsize=8)
    axs[2].grid(alpha=0.3)
    fig.suptitle('图 8-10（自绘仿真）习题 1：a=4、N=128、L=4 m 的啁啾信号 DFT'
                 '——吸收边界预乘对混叠的抑制', fontsize=12)
    fig.subplots_adjust(top=0.84, bottom=0.12, left=0.06, right=0.98,
                        wspace=0.25)
    fp = os.path.join(outdir, 'out_fig10_prob1.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    return e_raw, e_win


def prob2():
    """习题 2：式 (8.20)->(8.22) 的代数填充（推导打印）。"""
    print('===== 习题 2：式 (8.20) -> (8.22) 的代数填充（约束 3 与分段数无关） =====')
    print('    起点（式 8.20）：(1+Dz1/R)*d1 - w*Dz1/D1 <= d2 <= (1+Dz1/R)*d1 + w*Dz1/D1')
    print('    代入 d2 = (1-a2)*d1 + a2*dn、Dz1 = a2*Dz（式 8.21）：')
    print('      (1+a2*Dz/R)*d1 - w*a2*Dz/D1 <= (1-a2)*d1 + a2*dn')
    print('                                <= (1+a2*Dz/R)*d1 + w*a2*Dz/D1')
    print('    左：a2*dn >= a2*d1 + a2*Dz*d1/R - w*a2*Dz/D1  （移项合并 (1-a2)d1）')
    print('    右：a2*dn <= a2*d1 + a2*Dz*d1/R + w*a2*Dz/D1')
    print('    两边除以 a2（a2>0）——alpha_2 严格消去：')
    print('      (1+Dz/R)*d1 - w*Dz/D1 <= dn <= (1+Dz/R)*d1 + w*Dz/D1   （式 8.22）')
    print('    与式 (7.53) 完全相同 —— 约束 3 对任意分段数一致；')
    print('    数值复核见工程 1 验证 0（1000 组随机参数判定不一致 = 0 组）。')


def prob3(outdir):
    """习题 3：点源 100 km 的采样图，单次 vs 分段。"""
    print('===== 习题 3：点源 100 km -> 2 m 望远镜的采样图与单次/分段对比 =====')
    wvl = 1e-6
    Dz = 100e3
    D2 = 2.0
    # 第 7 章习题 2 的结论：delta1 = delta2 = 1 cm、N >= 1024（约束 4 主导；
    # 点源 D1->0 使约束 1 失效）
    d1 = d2 = 1e-2
    N = 1024
    print('    取 delta1 = delta2 = 1 cm、N = %d（第 7 章习题 2 的选点）' % N)
    print('    约束 2（D1->0）：N >= D2/(2*dn) + w*Dz/(2*d1*dn) = %.0f'
          % (D2 / (2 * d2) + wvl * Dz / (2 * d1 * d2)))
    print('    约束 4（单次）：N >= w*Dz/(d1*dn) = %.0f  -> 主导' % (wvl * Dz / (d1 * d2)))
    dzm = dz_max_partial(d1, d2, N, wvl)
    print('    分段视角：Delta_z_max = min^2*N/lambda = %.3e m >= Dz = 1e5 m'
          '（余量仅 %.0f%% —— N=1024 正是约束 4 逼出来的）-> 单次可行；'
          '分 2 段则每段约束 4 只需 N >= %.0f，余量 x4'
          % (dzm, (dzm / Dz - 1) * 100, wvl * Dz / 2 / (d1 * d2)))
    print('    对照书例（delta1=66.7um, delta_n=200um, N=128, Dz=2m）：'
          'Delta_z_max = 0.57 m << 2 m -> 必须分段')
    print('    结论：分段的收益由 min(delta1,deltan)^2*N/lambda 相对 Dz 的大小决定；'
          '点源长距场景 delta~cm 量级，约束 4 天然宽松')
    # 采样图：约束 2 的 log2 N 等高线（约束 1 对 D1->0 失效，不画）
    d1g = np.logspace(-3.3, -1.3, 300)
    dng = np.logspace(-3.3, -1.0, 300)
    D1g, Dng = np.meshgrid(d1g, dng)
    Nc2 = D2 / (2 * Dng) + wvl * Dz / (2 * D1g * Dng)      # 约束 2（几何项略去）
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    cs = ax.contourf(D1g, Dng, np.log2(np.maximum(Nc2, 1)),
                     levels=np.arange(6, 20, 2), cmap='Greys')
    cb = fig.colorbar(cs, ax=ax, shrink=0.9)
    cb.set_label('log2 N（约束 2，几何项略去）')
    # 约束 4 的单次 N=1024 界：w*Dz/(d1*dn) <= 1024 -> d1*dn >= w*Dz/1024
    prod_min = wvl * Dz / N
    d1line = np.logspace(-3.3, -1.3, 200)
    ax.plot(d1line, prod_min / d1line, 'b-', lw=2,
            label='约束 4（单次，N=1024）：$\\delta_1\\delta_n \\geq$ 1e-3 m$^2$')
    # 分段（例如 20 段，每段 5 km）后的约束 4 界
    ax.plot(d1line, prod_min / 20 / d1line, 'b--', lw=1.5,
            label='约束 4（分 20 段）：界线下移 20 倍')
    ax.plot(d1, d2, 'r*', ms=15, mec='k', zorder=6,
            label='选点 $\\delta_1=\\delta_2$=1 cm，N=1024')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plain_log(ax.xaxis)
    plain_log(ax.yaxis)
    ax.set_xlabel('$\\delta_1$ [m]')
    ax.set_ylabel('$\\delta_n$ [m]')
    ax.set_title('图 8-11（自绘）习题 3：点源 100 km 的采样图'
                 '（红* 上方约束 4 本来就不紧，分段无收益）', fontsize=10)
    ax.legend(fontsize=8, loc='lower left')
    ax.grid(alpha=0.3, which='both')
    fig.tight_layout()
    fp = os.path.join(outdir, 'out_fig11_prob3.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)


def annulus_ft(f, D_out, D_in):
    """环孔光瞳 P = circ(2r/D_out) - circ(2r/D_in) 的（径向）傅里叶变换。

    FT{circ(r/a)} = pi*a^2 * 2*J1(2*pi*a*f)/(2*pi*a*f)，代入 a = D/2：
        P~(f) = (pi/4)*[D_out^2*2J1(pi*D_out*f)/(pi*D_out*f)
                        - D_in^2*2J1(pi*D_in*f)/(pi*D_in*f)]
    f=0 处取极限 = (pi/4)(D_out^2 - D_in^2)。
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))

    def jinc_term(D, fv):
        u = np.pi * D * fv
        small = np.abs(u) < 1e-12                 # 注意 |u|：负 f 也要走正规分支
        res = np.empty_like(u)
        res[small] = D ** 2                       # u->0: 2J1(u)/u -> 1
        uu = u[~small]
        res[~small] = D ** 2 * 2 * j1(uu) / uu
        return res

    out = (np.pi / 4.0) * (jinc_term(D_out, f) - jinc_term(D_in, f))
    return out


def prob4(outdir):
    """习题 4：环形孔径聚焦束 100 km 的采样分析 + 仿真 vs 解析。"""
    print('===== 习题 4：环形孔径（D_out=1.5m, D_in=0.5m）聚焦 100 km =====')
    wvl = 1.3e-6
    Dz = 100e3
    D_out, D_in = 1.5, 0.5
    D1 = D_out
    # 源面：30 点跨外径 -> delta1 = 50 mm；观察面：聚焦斑 Airy 半径
    # 1.22*lambda*Dz/D_out = 1.22*1.3e-6*1e5/1.5 = 106 mm —— delta_n 细于它
    # 才能分辨主瓣与环结构，取 delta_n = 2.5 mm（约 40 点/Airy 半径，充裕）；
    # N 由约束 2 逼到 1024
    delta1 = D1 / 30
    deltan = 2.5e-3
    N = 1024
    D2 = 0.5                                             # ROI 直径 0.5 m
    # 约束 2：
    c2 = constraint2_N(delta1, deltan, D1, D2, wvl, Dz)
    print('    delta1 = %.1f mm（30 点跨外径）、delta_n = %.1f mm（Airy 半径 10.6 mm）、'
          'N = %d' % (delta1 * 1e3, deltan * 1e3, N))
    print('    约束 2（ROI 0.5m）：N >= %.1f = 2^%.2f -> N=%d 成立'
          % (c2, np.log2(c2), N))
    # 约束 3（聚焦源 R = -Dz）：
    lo3, hi3 = constraint3_d2(delta1, D1, wvl, Dz, -Dz)
    print('    约束 3（R=-Dz，聚焦到目标）：|1+Dz/R| = 0 < D2/D1 -> 免检'
          '（区间 [%.0f, %.0f] mm 恒含 delta_n）' % (lo3 * 1e3, hi3 * 1e3))
    # 约束 4 / 分段：
    dzm = dz_max_partial(delta1, deltan, N, wvl)
    nseg = int(np.ceil(Dz / dzm))
    print('    约束 4：Delta_z_max = min^2*N/lambda = %.0f m -> 至少 %d 段（%d 平面）'
          % (dzm, nseg, nseg + 1))
    print('    取 %d 段（每段 %.1f km，<= Delta_z_max）' % (nseg, Dz / nseg / 1e3))
    nseg_use = nseg
    z = (np.arange(1, nseg_use + 1)) * Dz / nseg_use
    # 源：环孔 x 聚焦相位 exp(-i k r1^2/(2 Dz))
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * delta1, n * delta1)
    r1 = np.sqrt(x1 ** 2 + y1 ** 2)
    ap = circle(x1, y1, D_out / 2) - circle(x1, y1, D_in / 2)
    src = ap * np.exp(-1j * np.pi / wvl / Dz * r1 ** 2)
    print('    仿真：%d 次分段传播（收缩网格 m=%.3f，网格跟随光束收缩）……'
          % (nseg_use, deltan / delta1))
    xn, yn, Uout = ang_spec_multi_prop_vac(src, wvl, delta1, deltan, z)
    # 解析：聚焦相位与传播核啁啾相消 -> 观察面 = e^{ikDz}e^{ikr^2/2Dz}/(i w Dz) x P~(f)
    sl = N // 2
    xn1 = xn[sl, :]
    f = xn1 / (wvl * Dz)
    Ptil = annulus_ft(f, D_out, D_in)
    U_an = np.exp(1j * np.pi / wvl / Dz * xn1 ** 2) / (1j * wvl * Dz) * Ptil
    I_num = np.abs(Uout[sl, :]) ** 2
    I_an = np.abs(U_an) ** 2
    pk = I_an.max()
    win = np.abs(xn1) <= D2 / 2
    eI = np.max(np.abs(I_num[win] - I_an[win])) / pk
    mph = win & (I_an > 0.05 * pk)
    eph = np.max(np.abs(np.angle(Uout[sl, :][mph] / U_an[mph])))
    print('    解析峰值辐照度 = %.3e；仿真/解析峰值比 = %.4f'
          % (pk, I_num.max() / pk))
    print('    ROI(+-0.25 m)内：|dI|/峰值 = %.2e；相位残差(|U|>5%%峰值) = %.2e rad'
          % (eI, eph))
    # 采样图（约束 2 等高线 + 选点）
    d1g = np.logspace(np.log10(delta1 / 5), np.log10(delta1 * 3), 260)
    dng = np.logspace(np.log10(deltan / 6), np.log10(deltan * 4), 260)
    D1g, Dng = np.meshgrid(d1g, dng)
    Nc2 = constraint2_N(D1g, Dng, D1, D2, wvl, Dz)
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.8))
    ax = axs[0]
    cs = ax.contourf(D1g * 1e3, Dng * 1e3, np.log2(np.maximum(Nc2, 1)),
                     levels=np.arange(4, 14), cmap='Greys')
    cb = fig.colorbar(cs, ax=ax, shrink=0.9)
    cb.set_label('log2 N（约束 2）')
    ax.plot(d1g * 1e3, constraint1_d2(d1g, D1, D2, wvl, Dz) * 1e3, 'k--', lw=1.6,
            label='约束 1')
    ax.plot(delta1 * 1e3, deltan * 1e3, 'r*', ms=14, mec='k', zorder=6,
            label='选点 (%.0f, %.1f) mm, N=%d' % (delta1 * 1e3, deltan * 1e3, N))
    ax.set_xscale('log')
    ax.set_yscale('log')
    plain_log(ax.xaxis)
    plain_log(ax.yaxis)
    ax.set_xlabel('$\\delta_1$ [mm]')
    ax.set_ylabel('$\\delta_n$ [mm]')
    ax.set_title('(a) 习题 4(a)：采样分析（D1=1.5m, D2=0.5m,\\n'
                 'Δz=100 km, λ=1.3μm；Δz$_{max}$=%.0f km -> %d 段）'
                 % (dzm / 1e3, nseg_use), fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which='both')
    ax = axs[1]
    ax.plot(xn1, I_an / pk, 'k-', lw=1.3, label='解析辐照度')
    ax.plot(xn1, I_num / pk, 'r.', ms=4, label='仿真（%d 段）' % nseg_use)
    ax.set_xlim(-0.1, 0.1)
    ax.set_xlabel('$x_n$ [m]')
    ax.set_ylabel('归一化辐照度')
    ax.set_title('(b) 习题 4(b)：$y_n$=0 切片辐照度（主瓣 ±0.1 m，'
                 'Airy 半径 1.22λΔz/D=106 mm）\n'
                 r'|dI|/pk=%.1e（ROI ±0.25 m）' % eI, fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.suptitle('图 8-12（自绘仿真）习题 4：环形孔径聚焦束传播 100 km'
                 '（λ=1.3μm，N=1024，%d 段分段传播）' % nseg_use, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fp = os.path.join(outdir, 'out_fig12_prob4.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)
    # 相位切片单独量化（解析相位含二次相位 + jinc 差的 π 翻转）；
    # 拟合取中央主瓣（|x| <= 8 mm，P~(f) 无符号翻转，unwrap 不滑档）
    central = np.abs(xn1) <= 8e-3
    ph_num = np.angle(Uout[sl, :])
    ph_an = np.angle(U_an)
    print('    相位切片：二次相位系数拟合（仿真）= %.4e /m^2，解析 = %.4e /m^2'
          % (np.polyfit(xn1[central] ** 2, np.unwrap(ph_num[central]), 1)[0],
             np.polyfit(xn1[central] ** 2, np.unwrap(ph_an[central]), 1)[0]))
    return eI, eph


def main(outdir=None):
    outdir = outdir or os.path.dirname(os.path.abspath(__file__))
    prob1(outdir)
    prob2()
    prob3(outdir)
    prob4(outdir)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
