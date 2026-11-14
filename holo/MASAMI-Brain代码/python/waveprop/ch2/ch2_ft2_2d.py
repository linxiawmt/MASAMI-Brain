# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 2 章 =====
二维 DFT：ft2 / ift2 与方形孔径的频谱（对应原书 2.6 节 Listings 2.5–2.6）

复用说明：
    base.py（python/fourier/common/）的 ft2_centered 与本脚本按原书定义的
    ft2 在偶数 N 下完全一致；sinc 也取自 base（归一化定义 sin(pi x)/(pi x)）。
    import 路径处理同 ch2_ft_basics.py。

示例：方形孔径 g(x,y) = rect(x/D) rect(y/D)，其连续 FT 为
    G(fx,fy) = D^2 sinc(D fx) sinc(D fy)
这正是第 4 章以后"夫琅禾费衍射图样 = 孔径 FT"的最小例子（方孔的远场
衍射图样）。相位在 sinc 过零点两侧翻转 0/pi，呈棋盘结构。

数值验证（运行自动打印）：
    (1) ift2(ft2(g)) 往返最大误差；
    (2) 二维 Parseval：sum|g|^2 dx^2 = sum|G|^2 dfx dfy；
    (3) 中心样本 G(0,0) = 孔径面积 D^2；
    (4) fy=0 剖面 vs 解析 D^2 sinc(D fx) 的最大相对偏差；
    (5) 原书 Listing 2.6 的 ift2 假定方阵（N = size(G,1) 同时用于两维），
        本脚本对非方阵验证"必须用 (M*dfx)*(N*dfy)"的修正版。
MATLAB/Python 差异：
    - MATLAB ift2 用 (N*delta_f)^2 且 N 取 size(G,1)——仅对方阵正确；
      numpy 侧本脚本写成对 M x N 通用的形式。
    - imshow 默认 origin='upper'（图像惯例），MATLAB imagesc + axis xy
      才与数学坐标一致；本脚本显式传 origin='lower'。
用法：python ch2_ft2_2d.py <输出目录>  保存 PNG；否则弹窗显示。
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
from base import sinc  # noqa: E402


# ===== 书 2.6 节 Listing 2.5 / 2.6：二维定标正确的 DFT / DIFT =====
def ft2(g, delta):
    """二维 DFT（原书 Listing 2.5）：G = fftshift(fft2(fftshift(g)))*delta^2"""
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g))) * delta ** 2


def ift2(G, delta_f, delta_f_y=None):
    """二维 DIFT（原书 Listing 2.6 的通用化）：
    g = ifftshift(ifft2(ifftshift(G))) * (M*df_y)*(N*df_x)。
    delta_f 是 x 方向频率间隔 1/(N*delta_x)；delta_f_y 是 y 方向的
    1/(M*delta_y)，方阵且两维等间距时与 delta_f 相同（默认）。
    原书按方阵写死 (N*delta_f)^2（MATLAB 的 N=size(G,1) 是行数），
    非方阵时行/列两维的 L 各不相同，必须分开乘。"""
    M, N = G.shape
    dfy = delta_f_y if delta_f_y is not None else delta_f
    return np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(G))) \
        * (M * dfy) * (N * delta_f)


def main(outdir=None):
    # ---- 网格与孔径 ----
    N = 128            # 每边采样点数（原书约定偶数 N；FFT 对 2 的幂最快）
    L = 4.0            # 空域边长 [m] -> delta = L/N = 31.25 mm
    D = 0.5            # 方形孔径边长 [m]（D/delta = 16，恰为整数个胞元）
    delta = L / N
    delta_f = 1.0 / (N * delta)
    n = np.arange(N) - N // 2
    # 胞元中心采样（第 1 章笔记的约定）：x = (n + 1/2)*delta，孔径边界
    # 恰好落在胞元边界上 -> 离散孔径与连续 rect 逐点相等，面积严格 = D^2。
    # 若用原书 Listing 2.3 的节点采样 x = n*delta，边界样本 |x| = D/2 会被
    # 整格计入（17 格），G(0,0) 多出 2*D*delta 量级的偏差。
    x = (n + 0.5) * delta
    f = n / (N * delta)
    X, Y = np.meshgrid(x, x)          # indexing='xy'：与 MATLAB meshgrid 同义
    FX, FY = np.meshgrid(f, f)

    g = ((np.abs(X) <= D / 2) & (np.abs(Y) <= D / 2)).astype(float)
    G = ft2(g, delta)

    # ---- 数值验证 ----
    print('===== 二维 DFT 验证（方孔 rect(x/D)rect(y/D)，D=%.2f m）=====' % D)
    g_rt = ift2(G, delta_f)
    print('(1) max|ift2(ft2(g)) - g| = %.2e' % np.max(np.abs(g_rt - g)))
    lhs = np.sum(np.abs(g) ** 2) * delta ** 2
    rhs = np.sum(np.abs(G) ** 2) * delta_f ** 2
    print('(2) Parseval: sum|g|^2*dx^2 = %.10f, sum|G|^2*df^2 = %.10f (偏差 %.2e)'
          % (lhs, rhs, abs(lhs - rhs) / lhs))
    i0 = N // 2
    print('(3) |G(0,0)| = %.6f，孔径面积 D^2 = %.6f（偏差 %.2e）'
          % (abs(G[i0, i0]), D ** 2, abs(abs(G[i0, i0]) - D ** 2)))
    prof = np.abs(G[i0, :])                      # fy = 0 剖面（取幅值：胞元中心
    prof_ana = np.abs(D ** 2 * sinc(D * f))      # 采样给谱附加 pi*f*delta 线性相位）
    ok = np.abs(f) <= 2 / D                      # 主瓣 + 前两个旁瓣内比较
    rel = np.max(np.abs(prof - prof_ana)[ok]) / D ** 2
    print('(4) fy=0 剖面（|fx| <= 2/D）vs 解析 D^2 sinc(D fx)：最大偏差 %.2e'
          '（相对峰值）' % rel)

    # 非方阵检查：M x N 矩形阵列（两维 delta 相同、频率间隔不同）上，
    # 原书 ift2 的 (size(G,1)*df)^2 定标错；正确做法是两维分别乘 M*df_y、N*df_x
    M2 = 64
    gx = np.exp(-np.pi * (x / 1.0) ** 2)         # x 向（N 点）
    gy = np.exp(-np.pi * (x[:M2] / 1.0) ** 2)    # y 向（M2 点）
    g2 = np.outer(gy, gx)                        # (M2, N) 各向异性高斯
    G2 = ft2(g2, delta)
    dfx = 1.0 / (N * delta)                      # x 向频率间隔
    dfy = 1.0 / (M2 * delta)                     # y 向频率间隔
    rt_book = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(G2))) \
        * (G2.shape[0] * delta_f) ** 2           # 原书写法：MATLAB size(G,1)=行数 M2
    rt_fix = ift2(G2, dfx, delta_f_y=dfy)
    print('(5) %dx%d 非方阵：原书 (行数*df_x)^2 往返误差 %.2e；'
          ' (M*df_y)(N*df_x) 往返误差 %.2e'
          % (M2, N, np.max(np.abs(rt_book - g2)), np.max(np.abs(rt_fix - g2))))

    # ---- 绘图 ----
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    im0 = axes[0].imshow(g, extent=[x[0], x[-1], x[0], x[-1]],
                         origin='lower', cmap='gray')
    axes[0].set_title('孔径 $g(x,y)=\\mathrm{rect}(x/D)\\mathrm{rect}(y/D)$')
    axes[0].set_xlabel('x [m]'); axes[0].set_ylabel('y [m]')
    fig.colorbar(im0, ax=axes[0], shrink=0.85)

    im1 = axes[1].imshow(20 * np.log10(np.abs(G) / np.abs(G).max() + 1e-12),
                         extent=[f[0], f[-1], f[0], f[-1]],
                         origin='lower', cmap='magma', vmin=-60, vmax=0)
    axes[1].set_title(r'$|G(f_x,f_y)|$ [dB]，$G=D^2\,\mathrm{sinc}\,\mathrm{sinc}$')
    axes[1].set_xlabel(r'$f_x$ [m$^{-1}$]'); axes[1].set_ylabel(r'$f_y$ [m$^{-1}$]')
    fig.colorbar(im1, ax=axes[1], shrink=0.85)

    axes[2].plot(f[ok], prof_ana[ok], color='0.6', lw=2, ls='--',
                 label=r'解析 $D^2\,\mathrm{sinc}(Df_x)$')
    axes[2].plot(f[ok], prof[ok], 'b.', ms=4, label='DFT 剖面 ($f_y$=0)')
    axes[2].set_xlabel(r'$f_x$ [m$^{-1}$]'); axes[2].set_ylabel('幅值')
    axes[2].set_title('剖面对比（最大偏差 %.1e）' % rel)
    axes[2].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.tight_layout()

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        fp = os.path.join(outdir, 'ft2_demo.png')
        fig.savefig(fp, dpi=150)
        print('已保存：%s' % fp)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
