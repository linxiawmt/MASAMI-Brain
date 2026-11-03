# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 4 章 =====
夫琅禾费衍射（原书 4.1 节，Listings 4.1-4.2 与图 4.1）

原书要点：
    夫琅禾费衍射积分（式 4.1）在 Delta_z > 2 D^2 / lambda（式 4.2，远场判据）时
    成立，它可以写成一次二维 FT（式 4.3），观察面坐标与源面空间频率按
        x2 = lambda*Delta_z*fx1,  y2 = lambda*Delta_z*fy1   （式 4.4）
    线性映射。数值实现 = 一次定标正确的 ft2（Listing 4.1 fraunhofer_prop）。

数值验证（运行本脚本全部自动打印）：
    (0) jinc 约定核对：原书解析式用 jinc(D*r2/(lambda*Dz))，其中
        jinc(u) = 2 J1(pi u)/(pi u)（u=0 处为 1，第一个零点在 u=1.22，
        即爱里环半径 1.22*lambda*Dz/D）——与 base.py 圆盘谱 jinc 的关系；
    (1) 书例复现（N=512, L=7.5 mm, D=1 mm, lambda=1 um, Dz=20 m）：
        轴上辐照度、传感器区域（±0.2 m）内数值 vs 解析的最大偏差；
    (2) 能量守恒（离散 Parseval）：sum|U2|^2*d2^2 == 源面孔径面积
        （入射单位振幅平面波），并对比解析场在观察窗内的真实能量
        —— 差值即被"卷回"窗内的混叠能量；
    (3) 采样不足的混叠表现：N 从 512 减到 64（d1 变粗 -> 观察窗
        x2max = lambda*Dz/(2*d1) 变窄），图样比窗口宽时边缘卷回；
    (4) 两条路线对照：直接按菲涅耳衍射积分（式 4.6）做二维数值求和
        （不用 FFT）与 fraunhofer_prop 的 FT 路线逐点对比。

MATLAB/Python 差异（对应 .m 文件内注释）：
    - MATLAB 数组 1 起、numpy 0 起：居中排列的零坐标样本位置不同，
      频率轴一律 fX = (arange(N) - N//2)/(N*d1) 生成（与原书
      (-N/2 : N/2-1)/(N*d1) 相同的升序居中排列）。
    - np.meshgrid 默认 indexing='xy'，与 MATLAB meshgrid 一致：
      x2 沿列变化、y2 沿行变化。
    - 原书 ft2 = fftshift(fft2(fftshift(g)))*delta^2（Listing 2.5）；
      base.py 的 ft2_centered 对输入先做 ifftshift——偶数 N 时
      fftshift 与 ifftshift 相同，两者逐元素一致。
用法：python ch4_fraunhofer.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib + scipy（J0/J1，能量与 jinc 核对用）；
      base.py（同仓库 python/fourier/common，circle/jinc/ft2_centered）。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0, j1

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ---- 复用《计算傅里叶光学》笔记的公共基础模块（同仓库） ----
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import circle, jinc, ft2_centered  # noqa: E402


# ===== 书 Listing 4.1：fraunhofer_prop =====
def fraunhofer_prop(Uin, wvl, d1, Dz):
    """夫琅禾费传播（原书 Listing 4.1）。

    Uout = exp(i k/(2 Dz) (x2^2+y2^2)) / (i lambda Dz) * ft2(Uin, d1)
    二次相位因子只留在积分外（积分内的已被远场近似扔掉）；
    e^{i k Dz} 轴上相位与原书一样不计算。
    返回 Uout（居中排列）、x2、y2（升序坐标，单位 m）。
    """
    N = Uin.shape[0]                       # assume square grid
    k = 2 * np.pi / wvl
    fX = (np.arange(N) - N // 2) / (N * d1)   # 源面频率轴（升序、居中排列）
    x2, y2 = np.meshgrid(wvl * Dz * fX, wvl * Dz * fX)  # 式 (4.4) 的线性映射
    Uout = np.exp(1j * k / (2 * Dz) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * Dz) * ft2_centered(Uin, d1)
    return Uout, x2, y2


def jinc_s(u):
    """原书解析式里的 jinc：jinc(u) = 2 J1(pi u)/(pi u)，jinc(0)=1。"""
    u = np.asarray(u, dtype=float)
    out = np.ones_like(u)
    nz = u != 0
    out[nz] = 2.0 * j1(np.pi * u[nz]) / (np.pi * u[nz])
    return out


def circ(X, Y, D):
    """直径 D 的圆孔（原书 circ(x,y,D) 约定，|r| <= D/2 取 1）。"""
    return circle(X, Y, D / 2.0)


def main(outdir=None):
    # ---- 原书 Listing 4.2 的参数 ----
    N = 512          # number of grid points per side
    L = 7.5e-3       # total size of the grid [m]
    d1 = L / N       # source-plane grid spacing [m]
    D = 1e-3         # diameter of the aperture [m]
    wvl = 1e-6       # optical wavelength [m]
    k = 2 * np.pi / wvl
    Dz = 20.0        # propagation distance [m]

    area = np.pi * D ** 2 / 4.0          # 圆孔面积 = FT 的直流值
    x1, y1 = np.meshgrid((np.arange(N) - N // 2) * d1,
                         (np.arange(N) - N // 2) * d1)
    Uin = circ(x1, y1, D)

    print('===== 验证 0：jinc 约定核对 =====')
    rho = np.linspace(0, 400, 4001)      # rho = r2/(lambda*Dz)
    diff = np.max(np.abs(jinc(rho, R=D / 2) - area * jinc_s(D * rho)))
    print('    max|base.jinc(rho, R=D/2) - 面积*jinc_s(D*rho)| = %.2e'
          % diff)
    print('    （即 base.jinc(rho,R=D/2) = 面积 × 2J1(pi D rho)/(pi D rho)，'
          '爱里图样）')
    u_first_zero = 3.8317 / np.pi
    print('    jinc_s 第一个零点 u = %.4f -> 爱里环半径 1.22*lambda*Dz/D = %.4f m'
          % (u_first_zero, u_first_zero * wvl * Dz / D))

    print('===== 验证 1 + 2：书例复现（N=512）与能量守恒 =====')
    Uout, x2, y2 = fraunhofer_prop(Uin, wvl, d1, Dz)
    d2 = wvl * Dz / L                    # 观察面采样间隔 = lambda*Dz*delta_f1
    x2max = wvl * Dz / (2 * d1)          # 观察窗半宽 = lambda*Dz*f1_max
    print('    d1 = %.3e m, d2 = lambda*Dz/L = %.4f m, 观察窗 ±%.3f m'
          % (d1, d2, x2max))
    U_th = np.exp(1j * k / (2 * Dz) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * Dz) * area * jinc_s(D * np.hypot(x2, y2) / (wvl * Dz))
    I_num, I_ana = np.abs(Uout) ** 2, np.abs(U_th) ** 2
    peak = (area / (wvl * Dz)) ** 2
    print('    轴上辐照度：数值 %.6f mW/m^2 vs 解析 %.6f mW/m^2（高 %.2f%%）'
          % (I_num[N // 2, N // 2] * 1e3, peak * 1e3,
             100 * (I_num[N // 2, N // 2] / peak - 1)))
    # 两类误差分开看：
    #   标度误差 = 数值整体比解析高（节点采样把孔径边界样本整格计入，
    #              直流样本偏大 -> 各处等比例抬高，见《计算傅里叶光学》
    #              笔记 ch2 的"胞元中心采样"坑）；
    #   形状误差 = 各自归一化后的残差（这才反映混叠/相位误差）。
    shape = np.abs(I_num / I_num.max() - I_ana / I_ana.max())
    sensor = np.abs(x2) <= 0.2           # 原书"0.4 m 宽的靶板传感器"
    print('    传感器区域（±0.2 m）：标度误差 %.2e，归一化形状误差 %.2e'
          % (np.max(np.abs(I_num - I_ana)[sensor]) / peak, shape[sensor].max()))
    imax = np.unravel_index(shape.argmax(), shape.shape)
    print('    形状误差最大值 %.2e 在 x2=%.3f m（第一爱里零点 0.0244 m 附近的陡坡）'
          % (shape.max(), x2[imax]))
    print('    边缘（|x2|>0.6 m）形状误差只有 %.1e——但那里的图样本身只有峰值的'
          ' 1e-6 倍，相对局部信号仍显著（原书"大范围显示时边缘开始出偏差"）'
          % shape[np.abs(x2) > 0.6].max())

    P_num = np.sum(I_num) * d2 ** 2      # 离散 Parseval：应恒等于孔径面积
    print('    能量：sum|U2|^2*d2^2 = %.6e（孔径面积 %.6e，相对偏差 %.2e）'
          % (P_num, area, abs(P_num - area) / area))
    # 解析场在观察窗内的真实能量（精细网格上数值积分）
    nf = 1501
    xf = np.linspace(-x2max, x2max, nf)
    Xf, Yf = np.meshgrid(xf, xf)
    I_fine = ((area / (wvl * Dz)) ** 2) * \
        jinc_s(D * np.hypot(Xf, Yf) / (wvl * Dz)) ** 2
    P_ana_win = np.trapezoid(np.trapezoid(I_fine, xf), xf)
    print('    解析场窗内能量 = %.6e（占总能量 %.2f%%）-> 混叠卷回约 %.2f%%'
          % (P_ana_win, 100 * P_ana_win / area, 100 * (1 - P_ana_win / area)))

    print('===== 验证 3：采样不足的混叠（N=64，观察窗变窄） =====')
    N2 = 64
    d1b = L / N2
    x1b, y1b = np.meshgrid((np.arange(N2) - N2 // 2) * d1b,
                           (np.arange(N2) - N2 // 2) * d1b)
    Uout_b, x2b, y2b = fraunhofer_prop(circ(x1b, y1b, D), wvl, d1b, Dz)
    d2b = wvl * Dz / L
    U_th_b = np.exp(1j * k / (2 * Dz) * (x2b ** 2 + y2b ** 2)) \
        / (1j * wvl * Dz) * area * \
        jinc_s(D * np.hypot(x2b, y2b) / (wvl * Dz))
    I_num_b, I_ana_b = np.abs(Uout_b) ** 2, np.abs(U_th_b) ** 2
    shape_b = np.abs(I_num_b / I_num_b.max() - I_ana_b / I_ana_b.max())
    print('    观察窗 ±%.4f m（爱里环半径 0.0244 m -> 窗内只有 ~3.5 个环）'
          % (wvl * Dz / (2 * d1b)))
    print('    归一化形状误差 max = %.3f（外环整个错位；标度误差另计 %.2f%%）'
          % (shape_b.max(), 100 * (I_num_b.max() / I_ana_b.max() - 1)))
    P_num_b = np.sum(I_num_b) * d2b ** 2
    xf_b = np.linspace(-wvl * Dz / (2 * d1b), wvl * Dz / (2 * d1b), 2001)
    Xfb, Yfb = np.meshgrid(xf_b, xf_b)
    I_fb = ((area / (wvl * Dz)) ** 2) * \
        jinc_s(D * np.hypot(Xfb, Yfb) / (wvl * Dz)) ** 2
    P_ana_b = np.trapezoid(np.trapezoid(I_fb, xf_b), xf_b)
    print('    能量账本：数值 sum|U2|^2*d2^2 = %.6e（仍=孔径面积！混叠不丢能量，'
          '只搬能量）' % P_num_b)
    print('              解析窗内真能量 = %.6e（占 %.2f%%）-> 卷回 %.2f%%'
          % (P_ana_b, 100 * P_ana_b / area, 100 * (1 - P_ana_b / area)))

    print('===== 验证 4：两条路线对照（直接积分 vs FT） =====')
    # 直接按菲涅耳衍射积分（式 4.6）做二维求和（y2=0 切片），不用 FFT：
    # U(x2,0) = 1/(i lambda Dz) * sum U1 * exp(i k/2Dz (x1^2+y1^2))
    #           * exp(-i k/Dz x1 x2) * d1^2
    # 只对孔径内的样本求和（孔径外 U1=0）。
    inside = Uin != 0
    xa, ya = x1[inside], y1[inside]
    quad = np.exp(1j * k / (2 * Dz) * (xa ** 2 + ya ** 2))   # 被远场近似扔掉的相位
    xs = x2[N // 2, :]                       # y2=0 行的 x2（升序）
    integ = np.einsum('i,ji->j', quad,
                      np.exp(-1j * k / Dz * xs[:, None] * xa[None, :]))
    U_dir = integ * d1 ** 2 / (1j * wvl * Dz)
    U_dir *= np.exp(1j * k / (2 * Dz) * xs ** 2)
    sl = N // 2
    diff_direct = np.max(np.abs(np.abs(U_dir) - np.abs(Uout[sl, :]))) \
        / np.max(np.abs(Uout[sl, :]))
    print('    y2=0 切片：max||U_direct|-|U_ft||/max|U_ft| = %.2e'
          % diff_direct)
    print('    （远场近似丢弃的二次相位在孔径边缘只有 %.3f rad，'
          '所以两条路线几乎重合）'
          % (k / (2 * Dz) * (D / 2) ** 2))

    if outdir:
        # ---- 图 fig-4-5：书例复现（原书图 4.1 的自绘版 + 二维图样） ----
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
        half = int(0.1 / d2) + 1
        sub = I_num[N // 2 - half:N // 2 + half, N // 2 - half:N // 2 + half]
        im = axes[0].imshow(sub,
                            extent=[-half * d2, half * d2, -half * d2, half * d2],
                            cmap='inferno', origin='lower')
        axes[0].set_title('(a) 观察面辐照度（±0.1 m，线性标度）')
        axes[0].set_xlabel(r'$x_2$ [m]')
        axes[0].set_ylabel(r'$y_2$ [m]')
        fig.colorbar(im, ax=axes[0], shrink=0.85)
        sl = slice(N // 2 - half, N // 2 + half)
        axes[1].plot(x2[N // 2, sl], I_ana[N // 2, sl],
                     color='0.6', lw=1.8, ls='--', label='解析')
        axes[1].plot(x2[N // 2, sl], I_num[N // 2, sl], 'b.', ms=4,
                     label='数值')
        axes[1].set_xlabel(r'$x_2$ [m]')
        axes[1].set_ylabel(r'辐照度 [W/m$^2$]')
        axes[1].set_title('(b) $y_2=0$ 切片：数值 vs 解析（原书图 4.1 复现）')
        axes[1].legend(fontsize=9)
        axes[2].plot(x2[N // 2, sl], shape[N // 2, sl], 'k.')
        axes[2].set_xlabel(r'$x_2$ [m]')
        axes[2].set_ylabel(r'归一化形状误差')
        axes[2].set_title('(c) 残差在爱里环陡坡处最大（直流标度差已扣除）')
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.tight_layout()
        fp = os.path.join(outdir, 'fig-4-5.png')
        fig.savefig(fp, dpi=150)
        print('\n已保存：%s' % fp)

        # ---- 图 fig-4-6：采样不足的混叠 ----
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
        hb = N2 // 2
        axes[0].plot(x2b[hb, :], I_ana_b[hb, :], color='0.6', lw=1.8,
                     ls='--', label='解析')
        axes[0].plot(x2b[hb, :], I_num_b[hb, :], 'b.', ms=4, label='数值（N=64）')
        axes[0].axvline(wvl * Dz / (2 * d1b), color='r', ls=':', lw=1)
        axes[0].axvline(-wvl * Dz / (2 * d1b), color='r', ls=':', lw=1)
        axes[0].set_xlabel(r'$x_2$ [m]')
        axes[0].set_ylabel(r'辐照度 [W/m$^2$]')
        axes[0].set_title('(a) N=64：观察窗 ±0.085 m，图样被卷回窗内')
        axes[0].legend(fontsize=9)
        half5 = int(0.1 / d2) + 1
        sl5 = slice(N // 2 - half5, N // 2 + half5)
        axes[1].semilogy(np.abs(x2[N // 2, sl5]), shape[N // 2, sl5] + 1e-16,
                         'k.', ms=3, label='N=512（窗 ±0.68 m）')
        xb_c = np.abs(x2b[hb, :])
        axes[1].semilogy(xb_c, shape_b[hb, :] + 1e-16, 'r.', ms=4,
                         label='N=64（窗 ±0.085 m）')
        axes[1].set_xlabel(r'$|x_2|$ [m]')
        axes[1].set_ylabel(r'归一化形状误差')
        axes[1].set_title('(b) 同一位置的误差：N=64（红）比 N=512（黑）大两个量级')
        axes[1].legend(fontsize=9)
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.tight_layout()
        fp = os.path.join(outdir, 'fig-4-6.png')
        fig.savefig(fp, dpi=150)
        print('已保存：%s' % fp)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
