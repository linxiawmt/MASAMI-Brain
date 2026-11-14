# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 2 章 =====
三个变换案例研究：sinc / 高斯 / 带二次相位的高斯（对应原书 2.5 节）

三种信号的带宽类型不同，采样间隔 delta 的确定方法也不同：
    1) sinc(a x)：硬限带，G = (1/a) rect(f/a)，带宽 |f| <= a/2，
       奈奎斯特判据 delta <= 1/a（恰好两倍余量时取 delta = 0.85/a，原书图 2.10）；
    2) exp(-pi (a x)^2)：软限带，G = (1/|a|) exp(-pi (f/a)^2)。
       按 1/e^2 频率取 delta_e2 = (1/(2a)) sqrt(pi/2)（原书式 2.25），
       频带边缘可见混叠；改按 p=0.01 电平取 delta（原书式 2.26/2.27）后混叠消失；
    3) exp(-pi (a x)^2) * exp(i pi (b x)^2)：二次相位把带宽拉宽，
       p 电平带宽由原书式 (2.31)：f_p = sqrt(-(a^2 + b^4/a^2)/pi * ln p)。
原书勘误（本脚本数值复现）：
    原书 2.5.3 节称 a=0.25、b=0.57 时 "f_{x,0.01} = 1.96 m^-1"，但按式 (2.31)
    代入 p=0.01 得 1.60 m^-1；1.96 m^-1 实为 p=0.001 的结果。
数值验证（运行自动打印）：
    - 案例 1：幅度波纹、频带边缘的相位 pi 跳变计数；
    - 案例 2：delta_e2 与 delta_0.01 的最左样本混叠抬升对比；
    - 案例 3：式 (2.31) 的 p 电平带宽 vs 数值谱、原书两个数值的归属，
      以及式 (2.29) 解析 FT 与直接积分的核对、图 2.14 三条 b 曲线的带宽。
MATLAB 注意：numpy 的 np.sinc(x) 是 sin(pi x)/(pi x)（归一化），
与 MATLAB Signal Processing Toolbox 的 sinc 同定义；为不添依赖，
本脚本自带 mysinc（与 base.py 的 sinc 相同，此处从 base 导入）。
用法：python ch2_case_studies.py  （纯数值输出，无图形）
依赖：numpy；base.py（同仓库 python/fourier/common）。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import sinc, coords  # noqa: E402


def ft1(g, delta):
    """居中一维 DFT（原书 Listing 2.1）"""
    return np.fft.fftshift(np.fft.fft(np.fft.fftshift(g))) * delta


def ft_ana_sinc(f, a):
    """sinc(a x) 的解析 FT：(1/a) rect(f/a)（原书式 2.22）"""
    return np.where(np.abs(f) < a / 2, 1.0 / a,
                    np.where(np.abs(f) == a / 2, 0.5 / a, 0.0))


def ft_ana_gauss(f, a):
    """exp(-pi (a x)^2) 的解析 FT（原书式 2.24）"""
    return 1.0 / abs(a) * np.exp(-np.pi * (f / a) ** 2)


def ft_ana_gauss_quad(f, a, b):
    """带二次相位高斯的解析 FT（原书式 2.29）"""
    c = a ** 2 - 1j * b ** 2
    return 1.0 / np.sqrt(c) * np.exp(-np.pi * f ** 2 / c)


def bandwidth_p(a, b, p):
    """p 电平带宽（原书式 2.31；b=0 时退化为式 2.27）"""
    return np.sqrt(-(a ** 2 + b ** 4 / a ** 2) / np.pi * np.log(p))


def main():
    print('===== 案例 1：sinc 信号（原书图 2.10：a=1.1, delta=0.85/a, N=32）=====')
    a = 1.1
    delta = 0.85 / a
    N = 32
    x, f = coords(N * delta, N)
    g = sinc(a * x)
    G = ft1(g, delta)
    G_a = ft_ana_sinc(f, a)
    band = np.abs(f) <= a / 2
    rip = np.max(np.abs(np.abs(G) - np.abs(G_a))[band]) * a
    edge = ~band
    ph_edge = np.angle(G[edge])
    frac_pi = np.mean(np.abs(np.abs(ph_edge) - np.pi) < 0.1)
    print('    带内幅度波纹 = %.2e（相对峰值 1/a=%.2f）' % (rip, 1 / a))
    print('    带外样本中相位≈pi 的比例 = %.0f%%（谱值应为 0，实为小负数）'
          % (100 * frac_pi))

    print('===== 案例 2：高斯信号（原书图 2.11/2.12：a=1）=====')
    a = 1.0
    f_e2 = a * np.sqrt(2 / np.pi)
    delta_e2 = 1.0 / (2 * f_e2)
    delta_01 = 1.0 / (2 * bandwidth_p(a, 0.0, 0.01))
    for tag, d in [('delta_e2  = %.4f' % delta_e2, delta_e2),
                   ('delta_0.01 = %.4f' % delta_01, delta_01)]:
        N = 32
        x, f = coords(N * d, N)
        G = ft1(np.exp(-np.pi * a * x ** 2), d)
        i_left = 0                      # 居中排列的最左样本 = 频带负边缘
        lift = G[i_left].real - ft_ana_gauss(f[i_left], a)
        print('    %s：最左样本 f=%.3f，混叠抬升 = %.4f（解析值 %.4f，抬升 %.0f%%）'
              % (tag, f[i_left], G[i_left].real,
                 ft_ana_gauss(f[i_left], a),
                 100 * lift / ft_ana_gauss(0, a)))

    print('===== 案例 3：二次相位高斯（原书图 2.13–2.15）=====')
    a, b = 0.25, 0.57
    # --- 式 (2.31) 的 p 电平带宽 vs 数值谱的 p 电平位置 ---
    for p in (0.01, 0.001):
        fp = bandwidth_p(a, b, p)
        # 数值：细网格上解析谱幅值降到 p 倍峰值的频率（验证式 2.29/2.31）
        f_fine = np.linspace(0, 4, 400001)
        A = np.abs(ft_ana_gauss_quad(f_fine, a, b))
        fp_num = np.interp(p, (A / A[0])[::-1], f_fine[::-1])
        print('    p=%-6g：式(2.31) f_p = %.4f m^-1，数值谱电平位置 = %.4f m^-1'
              % (p, fp, fp_num))
    print('    -> 原书正文"f_{x,0.01} = 1.96 m^-1"实为 p=0.001 的值（勘误）；'
          'p=0.01 应为 %.2f m^-1' % bandwidth_p(a, b, 0.01))
    # --- 式 (2.29) 解析 FT vs 直接数值积分 ---
    f0 = 0.9
    xx = np.linspace(-40, 40, 1_200_001)
    g_num = np.exp(-np.pi * (a * xx) ** 2) * np.exp(1j * np.pi * (b * xx) ** 2)
    G_num = np.trapezoid(g_num * np.exp(-2j * np.pi * f0 * xx), xx)
    rel = abs(G_num - ft_ana_gauss_quad(f0, a, b)) / abs(ft_ana_gauss_quad(f0, a, b))
    print('    式(2.29) 在 f=0.9 处 vs 直接积分：相对偏差 = %.2e' % rel)
    # --- 原书采样参数下的 DFT 精度 ---
    for tag, d, N in [('原书 delta=0.25 m（=p=0.001 判据）', 0.25, 40),
                      ('p=0.01 判据 delta=%.3f m' % (1 / (2 * bandwidth_p(a, b, 0.01))),
                       1 / (2 * bandwidth_p(a, b, 0.01)), 32)]:
        x, f = coords(N * d, N)
        g = np.exp(-np.pi * (a * x) ** 2) * np.exp(1j * np.pi * (b * x) ** 2)
        G = ft1(g, d)
        G_a = ft_ana_gauss_quad(f, a, b)
        ok = np.abs(G_a) > 1e-3 * np.abs(G_a).max()
        amp_err = np.max(np.abs(np.abs(G) - np.abs(G_a))[ok]) / np.abs(G_a).max()
        ph_err = np.max(np.abs(np.angle(G / G_a))[ok])
        print('    %s N=%d：带内幅度误差 %.2e，相位误差 %.2e rad'
              % (tag, N, amp_err, ph_err))
    # --- 原书图 2.14：a=0.33，b 增大拉宽带宽 ---
    print('    图 2.14 的带宽（a=0.33，p=0.01）：')
    for bb in (0.75, 1.5, 2.5):
        print('      b=%-4.2f -> f_p = %6.2f m^-1' % (bb, bandwidth_p(0.33, bb, 0.01)))


if __name__ == '__main__':
    main()
