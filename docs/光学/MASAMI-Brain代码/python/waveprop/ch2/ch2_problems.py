# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 2 章 =====
原书 2.7 节习题 1–5：五个信号在 a=1 与 a=10 两档下的 DFT vs 解析 FT

题目：对下列信号做 DFT，画出与解析 FT 的对比（a=1、a=10 各一次）：
    1) sinc(a x)              -> (1/a) rect(f/a)            [硬限带]
    2) exp(-pi (a x)^2)       -> (1/|a|) exp(-pi (f/a)^2)   [软限带]
    3) exp(-pi (a x)^2 + i pi (b x)^2), b=2
                              -> 1/sqrt(a^2-i b^2) exp(-pi f^2/(a^2-i b^2))
                                 [软限带，b 拉宽带宽]
    4) tri(a x)               -> (1/a) sinc^2(f/a)          [非限带，空域有限支撑]
    5) exp(-a |x|)            -> 2a/(a^2 + (2 pi f)^2)      [非限带，空域无限支撑]
采样策略（原书 2.5 节的案例研究就是这套方法的模板）：
    - sinc：奈奎斯特判据 delta = 0.85/a（留 15% 余量）；
    - 高斯类：p=0.01 电平带宽 -> delta = 1/(2 f_p)；
    - tri/exp(-a|x|)：谱按 1/f^2 衰减，取 |G| 降到峰值 p=0.01 的频率为
      f_max（tri 的 sinc^2 尾 -> f_max = a/(pi*sqrt(p))，洛伦兹尾 ->
      f_max = a sqrt(1/p-1)/(2 pi)），空域范围相应取振幅 < 0.01 处。
每个工况打印：delta、N、f_max 与 DFT 相对解析谱的最大偏差。
MATLAB 注意：tri 用 base 同款定义 |x|<=1 取 1-|x|；MATLAB 无内置 tri，
需自己写（见配套 .m 文件末尾的本地函数）。
用法：python ch2_problems.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib；base.py（同仓库 python/fourier/common）。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import sinc, tri, coords  # noqa: E402


def ft1(g, delta):
    """居中一维 DFT（原书 Listing 2.1）"""
    return np.fft.fftshift(np.fft.fft(np.fft.fftshift(g))) * delta


# ---- 五个信号与解析 FT ----
def f_sinc(x, a):
    return sinc(a * x)


def G_sinc(f, a):
    return np.where(np.abs(f) < a / 2, 1.0 / a,
                    np.where(np.abs(np.asarray(f)) == a / 2, 0.5 / a, 0.0))


def f_gauss(x, a):
    return np.exp(-np.pi * (a * x) ** 2)


def G_gauss(f, a):
    return 1.0 / abs(a) * np.exp(-np.pi * (f / a) ** 2)


def f_gauss_quad(x, a, b=2.0):
    return np.exp(-np.pi * (a * x) ** 2) * np.exp(1j * np.pi * (b * x) ** 2)


def G_gauss_quad(f, a, b=2.0):
    c = a ** 2 - 1j * b ** 2
    return 1.0 / np.sqrt(c) * np.exp(-np.pi * f ** 2 / c)


def f_expdec(x, a):
    return np.exp(-a * np.abs(x))


def G_expdec(f, a):
    return 2 * a / (a ** 2 + (2 * np.pi * f) ** 2)


def bandwidth_p_general(G_func, a, p, f_hi=4000.0):
    """数值求 |G| 降到峰值 p 倍的频率（用于无闭式带宽的信号）"""
    f = np.linspace(1e-9, f_hi, 2_000_001)
    A = np.abs(G_func(f, a)) / np.abs(G_func(0.0, a))
    return np.interp(p, A[::-1], f[::-1])


def main(outdir=None):
    a_list = [1.0, 10.0]
    p = 0.01
    b3 = 2.0                              # 习题 3 指定 b = 2
    problems = [
        ('1: sinc(ax)',            lambda x, a: f_sinc(x, a),      G_sinc),
        ('2: exp(-pi a^2 x^2)',    lambda x, a: f_gauss(x, a),     G_gauss),
        ('3: quad-phase, b=2',     lambda x, a: f_gauss_quad(x, a, b3),
                                   lambda f, a: G_gauss_quad(f, a, b3)),
        ('4: tri(ax)',             lambda x, a: tri(a * x),
                                   lambda f, a: sinc(f / a) ** 2 / a),
        ('5: exp(-a|x|)',          lambda x, a: f_expdec(x, a),    G_expdec),
    ]

    print('===== 习题 1–5：a=1 / a=10 的采样参数与 DFT 精度 =====')
    print('%-18s %-5s %-10s %-8s %-12s %-12s' %
          ('信号', 'a', 'delta', 'N', '带内误差', '带外误差'))
    rows = []
    for name, fx, Gx in problems:
        for a in a_list:
            if name.startswith('1'):
                delta = 0.85 / a                    # 奈奎斯特 + 15% 余量
            elif name.startswith('4'):
                f_max = a / (np.pi * np.sqrt(p))    # sinc^2 尾：sinc^2(u)=p
                delta = 1.0 / (2 * f_max)
            elif name.startswith('5'):
                f_max = bandwidth_p_general(G_expdec, a, p)
                delta = 1.0 / (2 * f_max)
            else:
                # 信号 2：f_p = sqrt(a^2/pi*ln(1/p))；信号 3：式 (2.31)（b=2）
                coeff = a ** 2 if name.startswith('2') \
                    else a ** 2 + b3 ** 4 / a ** 2
                f_max = np.sqrt(coeff / np.pi * np.log(1 / p))
                delta = 1.0 / (2 * f_max)
            # 空域范围：截到信号幅度 < 0.01（tri 本身有限支撑，取 1.2/a）
            if name.startswith('4'):
                half = 1.2 / a
            elif name.startswith('1'):
                half = 12 / a                       # sinc 尾按 1/x 截到 ~3%
            elif name.startswith('5'):
                half = np.log(100) / a              # exp(-a|x|) 幅度 0.01 处
            else:
                half = np.sqrt(np.log(100) / np.pi) / a   # 高斯幅度 0.01 处
            N = int(2 ** np.ceil(np.log2(2 * half / delta)))   # 2 的幂
            x, f = coords(N * delta, N)
            g = fx(x, a)
            G = ft1(g, delta)
            G_a = Gx(f, a)
            peak = np.max(np.abs(G_a))
            dev = np.abs(np.abs(G) - np.abs(G_a)) / peak
            if name.startswith('1'):     # 硬限带：带内看波纹，带外看 Gibbs 溢出
                inb = np.abs(f) < a / 2
                err_in = np.max(dev[inb])
                err_out = np.max(dev[~inb])
            else:
                err_in = np.max(dev)
                err_out = float('nan')
            rows.append((name, a, delta, N, x, f, np.abs(G), np.abs(G_a),
                         err_in, err_out))
            print('%-18s %-5.0f %-10.4g %-8d %-12.2e %-12.4g'
                  % (name, a, delta, N, err_in, err_out))

    # ---- 绘图：5 行 x 2 列（a=1 左，a=10 右）；显示范围 = 0.8/delta ≈ 1.6 倍带宽 ----
    fig, axes = plt.subplots(5, 2, figsize=(12, 16))
    for r, (name, _, _) in enumerate(problems):
        for c, a in enumerate(a_list):
            row = next(q for q in rows if q[0] == name and q[1] == a)
            _, _, delta, N, x, f, Ga, Gan, err_in, err_out = row
            ax = axes[r, c]
            m = np.abs(f) <= 0.8 / delta
            ax.plot(f[m], Gan[m], color='0.6', lw=1.8, ls='--', label='解析')
            ax.plot(f[m], Ga[m], 'b.', ms=3, label='DFT')
            err_txt = err_in if err_out != err_out else \
                '%.1e/%.1e' % (err_in, err_out)
            ax.set_title('%s   a=%g   (delta=%.3g, N=%d, 误差 %s)'
                         % (name, a, delta, N, err_txt), fontsize=10)
            if r == 4:
                ax.set_xlabel(r'$f_x$ [m$^{-1}$]')
            if r == 0 and c == 0:
                ax.legend(fontsize=9)
            ax.grid(alpha=0.25)
    for ax in axes[:, 0]:
        ax.set_ylabel('幅值')
    fig.tight_layout()

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        fp = os.path.join(outdir, 'problems.png')
        fig.savefig(fp, dpi=150)
        print('已保存：%s' % fp)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
