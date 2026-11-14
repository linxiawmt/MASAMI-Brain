# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 2 章 =====
离散化的三种失真与三条补救路线（对应原书 2.3、2.4 节与图 2.6–2.9）

连续 FT 对 g(x)=exp(-a|x|) ⇔ G(f)=2a/(a^2+(2*pi*f)^2)（原书式 (2.13)(2.14)，
a=10 m^-1，与原书图 2.6 参数一致）。DFT 链条的三步失真：
    1) 空域采样（乘 delta^-1 comb(x/delta)）  -> 频域周期复制 comb(delta*f)：混叠
    2) 空域截断（乘 rect(x/L)）              -> 频域卷积 L sinc(L f)：波纹+展宽
    3) 频域采样（乘 comb(L f)）              -> 空域虚拟周期延拓（周期 L）
其中第 1+2 步合并后的"连续频谱"有一个无需数值卷积的精确算法：
    S(f) = sum_{n=-N/2}^{N/2-1} g(n*delta) exp(-i 2 pi f n delta)
（这就是"采样 + 截断后的空域信号"的真实连续 FT；对它再在 f_m = m/(N*delta)
 上取样、乘 delta，就得到 DFT 输出——脚本内部对此做了自检。）

三条补救路线（原书图 2.7 / 2.8 / 2.9）：
    A 加大网格：delta 翻倍、N 不变     -> L=N*delta 翻倍，压波纹，但混叠变重
    B 加密采样：delta 减半、N 翻倍     -> L 不变，压混叠，但波纹不变
    C 双管齐下：delta 减半、L 也翻倍   -> N=4 倍，混叠与波纹同时下降
数值验证（运行自动打印）：
    - 混叠折返恒等式：6 m^-1 的余弦按 delta=1/8 m 采样与 2 m^-1 逐点相同；
    - delta*S(f_m) 与 ft(g) 的 DFT 输出逐 bin 一致（链条自洽）；
    - 四种策略的混叠/波纹/总误差对照表。
MATLAB 注意：本脚本无 MATLAB 专有函数；fftshift 语义与 MATLAB 相同
（偶数 N 时 fftshift == ifftshift）。
用法：python ch2_discretization.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


def g_func(x, a):
    """空域信号 exp(-a|x|)（原书式 2.13）"""
    return np.exp(-a * np.abs(x))


def G_func(f, a):
    """解析 FT：2a/(a^2+(2 pi f)^2)（原书式 2.14，峰值 2/a）"""
    return 2 * a / (a ** 2 + (2 * np.pi * f) ** 2)


def alias_only(f, a, delta, K=400):
    """只做第 1 步的频谱：G ⊗ comb(delta f) = (1/delta) Σ_k G(f + k/delta)。
    （原书式 2.15 的空域采样对 (1/delta)comb(x/delta) 在频域产生
    间隔 1/delta、幅度放大 1/delta 的周期复制；洛伦兹尾按 1/f^2 衰减，
    K=400 个周期后截断误差可忽略。）"""
    s = np.zeros_like(f, dtype=float)
    for k in range(-K, K + 1):
        s += G_func(f + k / delta, a)
    return s / delta


def chain_spectrum(f, g_samp, delta):
    """第 1+2 步的精确频谱 S(f) = Σ g(n delta) exp(-i2 pi f n delta)。
    g_samp 按居中排列（n = -N/2 .. N/2-1），返回复数谱（长度同 f）。"""
    n = np.arange(len(g_samp)) - len(g_samp) // 2
    terms = g_samp[np.newaxis, :] * np.exp(
        -2j * np.pi * f[:, np.newaxis] * (n * delta)[np.newaxis, :])
    return terms.sum(axis=1)


def ft1(g, delta):
    """居中一维 DFT（原书 Listing 2.1；等价于 base.py 的 ft1_centered）"""
    return np.fft.fftshift(np.fft.fft(np.fft.fftshift(g))) * delta


def main(outdir=None):
    a = 10.0
    # ---- 四种策略：(名称, delta, L, N) ----
    strategies = [
        ('基线（原书图 2.6）',        0.0375,  0.6, 16),
        ('A 加大网格：delta x2（图 2.7）', 0.075,  1.2, 16),
        ('B 加密采样：N x2（图 2.8）',   0.01875, 0.6, 32),
        ('C 双管齐下（图 2.9）',        0.01875, 1.2, 64),
    ]

    print('===== DFT 链条自检：delta*S(f_m) 与 DFT 输出逐 bin 对比 =====')
    f_fine = np.linspace(-20.0, 20.0, 8001)
    for name, delta, L, N in strategies:
        x = (np.arange(N) - N // 2) * delta
        g_samp = g_func(x, a)
        S = chain_spectrum(f_fine, g_samp, delta)
        fm = (np.arange(N) - N // 2) / (N * delta)      # DFT bin 频率
        S_bins = chain_spectrum(fm, g_samp, delta)
        G_dft = ft1(g_samp, delta)
        err = np.max(np.abs(delta * S_bins - G_dft))
        print('  %-24s max|delta*S - DFT| = %.2e' % (name, err))

    print('\n===== 混叠折返恒等式（原书图 2.5 第一组：f1=6, delta=1/8 m）=====')
    n_ali = np.arange(-16, 16)
    x_ali = n_ali / 8.0
    d_max = np.max(np.abs(np.cos(2 * np.pi * 6 * x_ali)
                          - np.cos(2 * np.pi * 2 * x_ali)))
    print('    max|cos(2*pi*6*x) 与 cos(2*pi*2*x) 的样本差| = %.1e（f_s=8, 折返 f2=2）'
          % d_max)

    print('\n===== 三种失真的量化（频带 |f| <= 1/(2 delta)，误差以解析峰值为基准）=====')
    print('  %-24s %-14s %-14s %-14s' %
          ('策略', '混叠误差', '波纹误差', 'DFT 总误差'))
    results = []
    for name, delta, L, N in strategies:
        x = (np.arange(N) - N // 2) * delta
        g_samp = g_func(x, a)
        fmax = 1.0 / (2 * delta)
        band = np.abs(f_fine) <= fmax
        G_true = G_func(f_fine, a)
        G_ali = alias_only(f_fine, a, delta)     # (1/delta)ΣG，放大标度
        S = chain_spectrum(f_fine, g_samp, delta).real   # 同为放大标度
        e_alias = np.max(np.abs(delta * G_ali - G_true)[band]) / G_func(0, a)
        e_rip = np.max(np.abs(delta * S - delta * G_ali)[band]) / G_func(0, a)
        fm = (np.arange(N) - N // 2) / (N * delta)
        e_dft = np.max(np.abs(ft1(g_samp, delta)
                              - G_func(fm, a))) / G_func(0, a)
        results.append((name, delta, L, N, f_fine, G_true, G_ali, S,
                        fm, ft1(g_samp, delta), e_alias, e_rip, e_dft))
        print('  %-24s %-14.2e %-14.2e %-14.2e'
              % (name, e_alias, e_rip, e_dft))

    # ---- 绘图：4 行 x 2 列（左列用原书图 (d)(f) 的 1/delta 放大标度）----
    fig, axes = plt.subplots(4, 2, figsize=(11.5, 13.5))
    for r, res in enumerate(results):
        name, delta, L, N, ff, Gt, Ga, S, fm, Gdft, *_ = res
        ax = axes[r, 0]
        ax.plot(ff, Gt / delta, color='0.65', lw=1.6, ls='--',
                label=r'解析 $G(f_x)/\delta$')
        ax.plot(ff, Ga, color='tab:blue', lw=1.0, ls=':',
                label=r'$G\otimes\mathrm{comb}(\delta f_x)$（混叠）')
        ax.plot(ff, S, color='k', lw=1.3,
                label=r'再 $\otimes L\,\mathrm{sinc}(Lf_x)$（+波纹）')
        ax.set_xlim(-20, 20)
        ax.set_ylim(-0.05, Ga.max() * 1.12)
        ax.set_title('%s   $\\delta$=%.4g m, L=%.2g m, N=%d'
                     % (name, delta, L, N), fontsize=10)
        ax.set_ylabel('幅值（$1/\\delta$ 放大标度）')
        if r == 0:
            ax.legend(fontsize=8, loc='upper right')
        if r == 3:
            ax.set_xlabel(r'$f_x$ [m$^{-1}$]')

        ax = axes[r, 1]
        markerline, stemlines, _ = ax.stem(fm, np.abs(Gdft))
        plt.setp(stemlines, lw=0.8)
        plt.setp(markerline, ms=3.5)
        ax.plot(ff, Gt, color='0.65', lw=1.6, ls='--', label='解析')
        ax.set_xlim(-20, 20)
        ax.set_ylim(-0.05, 0.32)
        ax.set_title('DFT 输出 vs 解析（总误差 %.1e）' % res[12], fontsize=10)
        if r == 0:
            ax.legend(fontsize=8, loc='upper right')
        if r == 3:
            ax.set_xlabel(r'$f_x$ [m$^{-1}$]')

    for ax in axes.ravel():
        ax.grid(alpha=0.25)
    fig.tight_layout()

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        fp = os.path.join(outdir, 'dft_remedies.png')
        fig.savefig(fp, dpi=150)
        print('\n已保存：%s' % fp)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
