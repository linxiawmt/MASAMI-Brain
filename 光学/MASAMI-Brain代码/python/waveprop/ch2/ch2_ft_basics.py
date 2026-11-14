# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 2 章 =====
DFT 基础：定标、重排、正逆变换配对（对应原书 2.1 节与 Listings 2.1–2.4）

复用说明：
    本脚本从 python/fourier/common/base.py 导入公共函数（sinc/gauss/coords）。
    base.py 属《计算傅里叶光学》(Voelz) 笔记的公共模块，其中的
    ft1_centered/ft2_centered 与本节按原书定义的 ft/ft2 在偶数 N 下完全等价
    （偶数 N 时 fftshift 与 ifftshift 相同）。import 路径处理：

        import os, sys
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "..", "fourier", "common"))
        from base import gauss, coords

数值验证（运行本脚本全部自动打印）：
    (1) 原书式 (2.6) 的 DFT 与"逐点黎曼和"一致（证明 fft 确实算的是那个和）；
    (2) 原书式 (2.5) 的重排公式按字面无法实现（原书笔误），修正公式与
        fftshift 逐元素一致；
    (3) 原书 Listing 2.3 的高斯 DFT 与解析 FT 的最大偏差；
    (4) 原书 Listing 2.4 平移高斯的相位斜率恢复 -2*pi*x0（与解析相位比对），
        并复现"网格不够大时相位只在中心附近正确"的原书论断；
    (5) ift(ft(g)) 往返 + Parseval 定理。
MATLAB/Python 差异：
    - MATLAB 数组 1 起、numpy 0 起：居中排列的"零坐标"样本在 MATLAB 的
      N/2+1 号，numpy 的 N//2 号。坐标轴一律用 x = (n - N/2)*delta 生成，
      与排列解耦。
    - numpy.fft.fft 正变换不除 N（归一化在 ifft），MATLAB 相同——
      定标因子 delta / (N*delta_f) 都要自己乘。
用法：直接运行打印验证；`python ch2_ft_basics.py` 无图形输出。
依赖：numpy（不需要 matplotlib）。
"""
import os
import sys

import numpy as np

# ---- 复用《计算傅里叶光学》笔记的公共基础模块（同仓库） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import gauss, coords  # noqa: E402


# ===== 书 2.1.3 节 Listing 2.1 / 2.2：定标正确的 DFT / DIFT =====
def ft(g, delta):
    """一维 DFT（原书 Listing 2.1）：
    G = fftshift(fft(fftshift(g))) * delta
    输入 g 按居中排列（零坐标样本在中间，偶数 N 时 numpy 索引 N//2），
    输出 G 同样居中排列；乘 delta 把求和升级为连续 FT 的近似。"""
    return np.fft.fftshift(np.fft.fft(np.fft.fftshift(g))) * delta


def ift(G, delta_f):
    """一维 DIFT（原书 Listing 2.2）：
    g = ifftshift(ifft(ifftshift(G))) * length(G) * delta_f
    逆变换的定标因子是 N*delta_f = 1/delta（保持 ft/ift 互逆）。"""
    return np.fft.ifftshift(np.fft.ifft(np.fft.ifftshift(G))) \
        * len(G) * delta_f


def main():
    print('===== 验证 1：式 (2.6) 的 DFT 就是重排后的黎曼和 =====')
    N, delta = 8, 0.5
    n = np.arange(N) - N // 2                      # n = -N/2 .. N/2-1
    x = n * delta
    g = gauss(x / 2.0)                             # 任意测试函数（居中排列）
    m = np.arange(N) - N // 2
    # 逐点按式 (2.4) 的指数计算（未重排的对称形式）
    G_sum = np.array([delta * np.sum(g * np.exp(-2j * np.pi * mi * n / N))
                      for mi in m])
    G_ft = ft(g, delta)
    print('    max|G_ft - G_黎曼和| = %.2e' % np.max(np.abs(G_ft - G_sum)))

    print('===== 验证 2：式 (2.5) 重排公式（原书笔误）=====')
    # 原书 (2.5) 第一支：g[(n'+N/2)*delta], n'=1..N/2+1 —— n'=1 时坐标
    # (1+N/2)*delta 超出采样范围 [-N/2, N/2-1]*delta，按字面无法取到样本；
    # 第二支 g[(n'-N-2)*delta] 整体错一个序号（末端只到 -2*delta）。
    # 正确的重排（与 fftshift 逐元素一致，也与原书图 2.1 一致）：
    #     g_tilde[n'] = g[(n'-1)*delta]        n' = 1 .. N/2   （1-based）
    #     g_tilde[n'] = g[(n'-1-N)*delta]      n' = N/2+1 .. N
    # 0-based numpy 写法：g_tilde = fftshift(g)（居中 -> 零坐标在最前）。
    book_idx = n[0]  # 仅占位避免未用警告
    del book_idx
    g_tilde_correct = np.concatenate([g[N // 2:], g[:N // 2]])
    g_tilde_fftshift = np.fft.fftshift(g)
    print('    修正公式 vs fftshift 逐元素最大差 = %.2e'
          % np.max(np.abs(g_tilde_correct - g_tilde_fftshift)))
    ok = ((n + N / 2) * delta <= x[-1] + 1e-12).mean()
    print('    原书第一支可取到的样本比例（N=%d）= %.0f%%（其余越界）'
          % (N, 100 * ok))

    print('===== 验证 3：Listing 2.3 高斯 DFT vs 解析 FT =====')
    # 原书参数：L = 5, N = 32, a = 1
    L, N, a = 5.0, 32, 1.0
    delta = L / N
    x, f = coords(L, N)          # x = n*delta（居中），f = n/(N*delta)
    g_samp = gauss(a * x)        # exp(-pi*(a*x)^2)
    G_dft = ft(g_samp, delta)
    G_ana = 1.0 / a * gauss(f / a)   # (1/|a|) exp(-pi*(f/a)^2)，式 (2.24)
    err = np.abs(G_dft - G_ana)
    i0 = np.argmin(np.abs(f))    # f = 0 的样本（偶数 N 时位于 N//2）
    print('    max|G_dft - G_ana| = %.2e，其中 f=0 处偏差 = %.2e'
          % (err.max(), err[i0]))

    print('===== 验证 4：Listing 2.4 平移高斯的相位 =====')
    # 原书参数：L = 10, N = 64, x0 = 5*delta, a = 1
    for (L, N) in [(5.0, 32), (10.0, 64)]:
        delta = L / N
        x0 = 5 * delta
        x, f = coords(L, N)
        g_samp = gauss(x - x0)                     # 平移后的高斯
        G_dft = ft(g_samp, delta)
        # 解析 FT：exp(-i 2 pi x0 f) * exp(-pi f^2)，用 -pi*f^2/a（a=1）
        G_ana = np.exp(-2j * np.pi * x0 * f) * gauss(f)
        # 只在幅度显著处（|G| > 1e-3 * 峰值）比较展开相位
        amp_ok = np.abs(G_dft) > 1e-3 * np.abs(G_dft).max()
        ph_err = np.angle(G_dft / G_ana)[amp_ok]
        print('    L=%-4.1f N=%d：幅度显著样本的相位误差 max = %.2e rad'
              % (L, N, np.max(np.abs(ph_err))))

    print('===== 验证 5：ft/ift 往返 + Parseval =====')
    L, N = 5.0, 64
    delta = L / N
    delta_f = 1.0 / (N * delta)
    x, f = coords(L, N)
    g = gauss(x / 1.5) * (1 + 0.3 * np.cos(2 * np.pi * 0.7 * x))
    g_rt = ift(ft(g, delta), delta_f)
    print('    max|ift(ft(g)) - g| = %.2e' % np.max(np.abs(g_rt - g)))
    lhs = np.sum(np.abs(g) ** 2) * delta
    rhs = np.sum(np.abs(ft(g, delta)) ** 2) * delta_f
    print('    Parseval: sum|g|^2*dx = %.10f, sum|G|^2*df = %.10f, 相对偏差 %.2e'
          % (lhs, rhs, abs(lhs - rhs) / lhs))


if __name__ == '__main__':
    main()
