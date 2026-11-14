# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 6 章（工程 2）=====
点源模型：sinc 模型与 sinc-Gaussian 模型（原书 6.6 节，Listings 6.8-6.9）

原书要点：
    真点源 = Dirac delta（式 6.82），谱无限宽 —— 离散网格必须用带限模型：
      Martin-Flatté（式 6.83）：exp(-r^2/2sigma^2) exp(-i r^2/2x0^2)；
      Flatté et al.（式 6.84）：exp(-r^2/2sigma^2) cos^2(r^2/2rho^2)；
      本书做法（逆向设计）：要求观测面给出"窗口加权的抛物面相位球面波"
      （式 6.85/6.86），反解源面模型（式 6.87-6.92）：
        方窗（宽 D，中心 r_c）：sinc 模型（式 6.92，Listing 6.8）；
        乘 Gaussian 抑制旁瓣 -> sinc-Gaussian 模型（Listing 6.9），
        观测面相位朝边缘的混叠明显减轻（原书图 6.7-6.12）。

原书勘误（本脚本按修正实现并数值验证，详见笔记）：
    式 (6.92) 的两个指数因子符号印反：e^{+i k r_c^2/(2*Dz)} e^{-i (k/Dz) r_c.r1}
    应为 e^{-i k r_c^2/(2*Dz)} e^{+i (k/Dz) r_c.r1}。
    修正后模型恰为居中模型的平移：U_pt(r1) = A (D/(lambda*Dz))^2
      * exp(-i k |r1-r_c|^2/(2 Dz)) * sinc[D(x1-xc)/(lambda Dz)] sinc[D(y1-yc)/(lambda Dz)]
    （Fresnel 核平移不变 -> 传播后给出窗口中心在 +r_c 的目标场；书印版本的
    额外相位斜坡等价于倾斜入射，把窗口推到 -r_c —— 验证 4 实测质心反号。）

数值验证（运行本脚本全部自动打印）：
    (0) 源面模型峰值辐照度：A^2 arg^4 = (lambda Dz)^2 (D/(lambda Dz))^4
        = D^4/(lambda Dz)^2 = 4.096 kW/m^2（对照原书图 6.7 的 ~4.1）；
    (1) Listing 6.8 复现：sinc 模型传播（ang_spec）-> 观测面辐照度
        窗内平顶 ~1 W/m^2（图 6.8）、相位 vs 解析抛物线（图 6.9）；
        平顶度与相位误差分窗内/窗外统计；
    (2) Listing 6.9 复现：sinc-Gaussian 模型 -> 同样统计（图 6.10-6.12），
        量化"高斯因子压旁瓣、减相位混叠"；
    (3) 能量守恒（Parseval）抽查；
    (4) 勘误 (6.92) 对照：离轴窗口 r_c = (1.2, 0.8) mm（δ1 的整数倍），
        书印模型 vs 修正模型传播后的辐照度质心（预言：书印 -r_c，修正 +r_c）、
        修正模型 vs 目标场（式 6.86）的复场偏差。

MATLAB/Python 差异（对应 .m 文件内注释）：
    - numpy 0-based / MATLAB 1-based：居中排列的零坐标样本位置不同；
    - cart2pol -> hypot/atan2；sinc 两语言均为归一化定义 sin(pi x)/(pi x)
      （MATLAB 需 Signal Processing Toolbox）。
用法：python ch6_point_source.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib；ch6_fresnel_prop.py（ang_spec_prop/ift2，同目录）；
      base.py（rect/ft2_centered，同仓库 python/fourier/common）。
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import rect  # noqa: E402
from ch6_fresnel_prop import ang_spec_prop, ift2, ft2  # noqa: E402  (ft2 未直接用，保持一致导入)


# ===== 书 Listing 6.8：sinc 模型点源（居中窗口，式 6.92 的 rc=0 情形） =====
def pt_source_sinc(x1, y1, D, wvl, Dz, A=None):
    """sinc 模型点源（式 6.92，rc = 0；Listing 6.8）。
    A 缺省取 lambda*Dz，使观测面窗内振幅为 1（原书"sets field amplitude
    to 1 in obs plane"）。"""
    if A is None:
        A = wvl * Dz
    arg = D / (wvl * Dz)
    r1sq = x1 ** 2 + y1 ** 2
    return A * arg ** 2 * np.exp(-1j * 2.0 * np.pi / wvl / (2.0 * Dz) * r1sq) \
        * np.sinc(arg * x1) * np.sinc(arg * y1)


# ===== 书 Listing 6.9：sinc-Gaussian 模型点源 =====
def pt_source_sinc_gauss(x1, y1, D, wvl, Dz, A=None):
    """sinc-Gaussian 模型点源（Listing 6.9：sinc 模型乘 exp(-(arg/4 r1)^2)）。"""
    if A is None:
        A = wvl * Dz
    arg = D / (wvl * Dz)
    r1sq = x1 ** 2 + y1 ** 2
    return A * arg ** 2 * np.exp(-1j * 2.0 * np.pi / wvl / (2.0 * Dz) * r1sq) \
        * np.sinc(arg * x1) * np.sinc(arg * y1) \
        * np.exp(-(arg / 4.0 * np.sqrt(r1sq)) ** 2)


# ===== 式 (6.92) 的离轴推广（书印 vs 修正） =====
def pt_source_offaxis(x1, y1, D, wvl, Dz, rc, corrected=True, A=None):
    """窗口中心 r_c = (xc, yc) 的方窗（宽 D）sinc 模型点源。

    corrected=True（本笔记修正版）：U_pt = A arg^2 e^{-i k |r1-rc|^2/(2 Dz)}
        * sinc[arg(x1-xc)] sinc[arg(y1-yc)] —— 恰为居中模型的平移；
    corrected=False（书印 6.92）：U_pt = A arg^2 e^{-i k r1^2/(2Dz)}
        * e^{+i k rc^2/(2Dz)} e^{-i (k/Dz) rc.r1} * sinc[...]（指数符号印反）。"""
    if A is None:
        A = wvl * Dz
    arg = D / (wvl * Dz)
    k = 2.0 * np.pi / wvl
    xc, yc = rc
    sinc_env = np.sinc(arg * (x1 - xc)) * np.sinc(arg * (y1 - yc))
    if corrected:
        ph = -1j * k * ((x1 - xc) ** 2 + (y1 - yc) ** 2) / (2.0 * Dz)
    else:
        ph = (-1j * k * (x1 ** 2 + y1 ** 2) / (2.0 * Dz)
              + 1j * k * (xc ** 2 + yc ** 2) / (2.0 * Dz)
              - 1j * k / Dz * (xc * x1 + yc * y1))
    return A * arg ** 2 * np.exp(ph) * sinc_env


def target_windowed_wave(x2, y2, D, wvl, Dz, rc, A=None):
    """目标观测面场（式 6.86/6.90，省略全局相位 e^{ik Dz}）：
    U = (A/(i lambda Dz)) W(r2 - rc) e^{i k |r2-rc|^2/(2 Dz)}，
    方窗 W = rect((x2-xc)/D) rect((y2-yc)/D)；A 缺省取 lambda*Dz
    （窗内振幅 1，与点源模型的振幅定标一致）。"""
    if A is None:
        A = wvl * Dz
    k = 2.0 * np.pi / wvl
    xc, yc = rc
    W = rect((x2 - xc) / D) * rect((y2 - yc) / D)
    return A * W * np.exp(1j * k * ((x2 - xc) ** 2 + (y2 - yc) ** 2)
                          / (2.0 * Dz)) / (1j * wvl * Dz)


def phase_err_stats(U, Uan, r2sq, r_in, r_out1, r_out2):
    """包裹相位差（去均值）的窗内/窗外 RMS 统计。"""
    dph = np.angle(U * np.conj(Uan))
    dph -= np.mean(dph[r2sq <= r_in ** 2])
    m_in = r2sq <= r_in ** 2
    m_out = (r2sq >= r_out1 ** 2) & (r2sq <= r_out2 ** 2)
    return (np.sqrt(np.mean(dph[m_in] ** 2)),
            np.sqrt(np.mean(dph[m_out] ** 2)))


def main(outdir=None):
    # ---- 原书 Listing 6.8/6.9 的参数 ----
    D = 8e-3            # 观测窗口直径 [m]
    wvl = 1e-6          # optical wavelength [m]
    Dz = 1.0            # propagation distance [m]
    arg = D / (wvl * Dz)
    d1 = 1.0 / (10 * arg)   # source-plane grid spacing（中心瓣内 10 个采样）
    d2 = D / 100            # observation-plane grid spacing
    N = 1024
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    x2, y2 = np.meshgrid(n * d2, n * d2)
    r2sq = x2 ** 2 + y2 ** 2

    print('===== 验证 0：源面模型峰值辐照度 =====')
    pt = pt_source_sinc(x1, y1, D, wvl, Dz)
    ptg = pt_source_sinc_gauss(x1, y1, D, wvl, Dz)
    I_pt = np.abs(pt) ** 2
    I_ptg = np.abs(ptg) ** 2
    print('    A = lambda*Dz = %.2e m；arg = D/(lambda Dz) = %.0f /m；'
          'delta1 = %.2f um（中心瓣 10 采样/半瓣），delta2 = %.0f um'
          % (wvl * Dz, arg, d1 * 1e6, d2 * 1e6))
    print('    sinc 模型：源面峰值辐照度 = %.3f kW/m^2（解析 A^2 arg^4 = %.3f，'
          '原书图 6.7 峰值 ~4.1）'
          % (I_pt.max() / 1e3, (wvl * Dz) ** 2 * arg ** 4 / 1e3))
    print('    sinc-Gauss 模型：源面峰值辐照度 = %.3f kW/m^2（高斯因子在 r=0 处为 1，'
          '压的是旁瓣不是峰值）' % (I_ptg.max() / 1e3))

    print('===== 验证 1/2：传播与观测面统计（原书图 6.8/6.9 与 6.11/6.12） =====')
    _, _, U_ob = ang_spec_prop(pt, wvl, d1, d2, Dz)
    _, _, U_obg = ang_spec_prop(ptg, wvl, d1, d2, Dz)
    I_ob, I_obg = np.abs(U_ob) ** 2, np.abs(U_obg) ** 2
    # 能量守恒
    E_in = np.sum(np.abs(pt) ** 2) * d1 ** 2
    E_out = np.sum(I_ob) * d2 ** 2
    print('    Parseval：sum|U|^2 d2^2 / sum|pt|^2 d1^2 = %.6f'
          % (E_out / E_in))
    # 平顶度：窗内（|r2| <= D/2，去 5% 边缘）
    m_win = r2sq <= (0.95 * D / 2) ** 2
    for tag, Iob in (('sinc', I_ob), ('sinc-Gauss', I_obg)):
        print('    %s：窗内辐照度 mean = %.4f W/m^2，std/mean = %.2e'
              % (tag, Iob[m_win].mean(), Iob[m_win].std() / Iob[m_win].mean()))
    # 相位误差：窗内 vs 窗外（混叠区）
    Uan = (np.exp(1j * 2.0 * np.pi / wvl / (2.0 * Dz) * r2sq)
           / (1j * wvl * Dz))          # 目标：W=1（窗内）的抛物线相位；窗外幅度无意义
    r_in, r_out1, r_out2 = 0.95 * D / 2, 1.10 * D / 2, 6e-3
    for tag, U in (('sinc', U_ob), ('sinc-Gauss', U_obg)):
        e_in, e_out = phase_err_stats(U, Uan, r2sq, r_in, r_out1, r_out2)
        m_ring = (r2sq >= (1.10 * D / 2) ** 2) & (r2sq <= (6e-3) ** 2)
        print('    %s：相位误差 RMS 窗内 = %.4f rad；窗外 D/2~6 mm 环带 = %.4f rad；'
              '环带峰值 |U| = %.3f（窗内幅值 1）'
              % (tag, e_in, e_out, np.abs(U[m_ring]).max()))

    print('===== 验证 3：勘误 (6.92)：离轴窗口 r_c = (1.2, 0.8) mm =====')
    # D = 4 mm 窗，r_c 取 delta1 的整数倍（离散平移无插值误差）
    Dw = 4e-3
    argw = Dw / (wvl * Dz)
    d1w = 1.0 / (10 * argw)
    rc = (48 * d1w, 32 * d1w)          # = (1.2 mm, 0.8 mm)
    x1w, y1w = np.meshgrid(n * d1w, n * d1w)
    x2w, y2w = np.meshgrid(n * d2, n * d2)
    pt_ok = pt_source_offaxis(x1w, y1w, Dw, wvl, Dz, rc, corrected=True)
    pt_bad = pt_source_offaxis(x1w, y1w, Dw, wvl, Dz, rc, corrected=False)
    _, _, U_ok = ang_spec_prop(pt_ok, wvl, d1w, d2, Dz)
    _, _, U_bad = ang_spec_prop(pt_bad, wvl, d1w, d2, Dz)
    I_ok, I_bad = np.abs(U_ok) ** 2, np.abs(U_bad) ** 2
    # 辐照度质心（只取显著区域，避免噪声主导）
    def centroid(I):
        w = I > 0.05 * I.max()
        return ((x2w[w] * I[w]).sum() / I[w].sum(),
                (y2w[w] * I[w]).sum() / I[w].sum())
    c_ok = centroid(I_ok)
    c_bad = centroid(I_bad)
    print('    目标：窗口中心 r_c = (%.2f, %.2f) mm' % (rc[0] * 1e3, rc[1] * 1e3))
    print('    修正模型：传播后辐照度质心 = (%.2f, %.2f) mm —— 落在 +r_c'
          % (c_ok[0] * 1e3, c_ok[1] * 1e3))
    print('    书印模型：传播后辐照度质心 = (%.2f, %.2f) mm —— 落在 -r_c'
          '（额外相位斜坡 = 倾斜入射，窗口被推到镜像位置）'
          % (c_bad[0] * 1e3, c_bad[1] * 1e3))
    # 修正模型 vs 目标场（式 6.86）复场对照（窗内 + 内部 80% 区）
    U_tgt = target_windowed_wave(x2w, y2w, Dw, wvl, Dz, rc)
    m_t = np.abs(U_tgt) > 0.5 * np.abs(U_tgt).max()
    m_in = (((x2w - rc[0]) / (Dw / 2)) ** 2 + ((y2w - rc[1]) / (Dw / 2)) ** 2) <= 0.8 ** 2
    for tag, msk in (('全窗', m_t), ('内部 80%', m_t & m_in)):
        ratio = U_ok[msk] / U_tgt[msk]
        print('    修正模型 vs 目标场（式 6.86，%s）：|比值-1| max = %.2e；'
              '相位残差(去均值) max = %.2e rad'
              % (tag, np.max(np.abs(np.abs(ratio) - 1.0)),
                 np.max(np.abs(np.angle(ratio) - np.mean(np.angle(ratio))))))

    if not outdir:
        plt.show()
        return

    # ---- 图 6-14（自绘）：点源模型总览 ----
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.6))
    sl = N // 2
    x1mm = x1[0, :] * 1e3
    m1 = np.abs(x1mm) <= 0.75
    axs[0, 0].plot(x1mm[m1], I_pt[sl, m1] / 1e3, 'k-', lw=1.2)
    axs[0, 0].set_title('(a) sinc 模型点源（源面辐照度切片，Listing 6.8）\n'
                        '峰值 %.2f kW/m$^2$，中心瓣半宽 %.3f mm'
                        % (I_pt.max() / 1e3, 1 / arg * 1e3))
    axs[0, 0].set_xlabel('x1 [mm]'); axs[0, 0].set_ylabel('辐照度 [kW/m$^2$]')
    axs[0, 0].grid(alpha=0.3)
    axs[0, 1].plot(x1mm[m1], I_ptg[sl, m1] / 1e3, 'k-', lw=1.2)
    axs[0, 1].set_title('(b) sinc-Gaussian 模型（Listing 6.9）\n'
                        '峰值 %.2f kW/m$^2$（高斯因子压旁瓣）'
                        % (I_ptg.max() / 1e3))
    axs[0, 1].set_xlabel('x1 [mm]'); axs[0, 1].set_ylabel('辐照度 [kW/m$^2$]')
    axs[0, 1].grid(alpha=0.3)
    x2mm = x2[0, :] * 1e3
    m2 = np.abs(x2mm) <= 8
    axs[0, 2].plot(x2mm[m2], I_ob[sl, m2], 'b-', lw=1.0,
                   label='sinc 模型')
    axs[0, 2].plot(x2mm[m2], I_obg[sl, m2], 'r-.', lw=1.0,
                   label='sinc-Gauss 模型')
    axs[0, 2].axvline(D / 2 * 1e3, color='0.5', ls=':')
    axs[0, 2].axvline(-D / 2 * 1e3, color='0.5', ls=':')
    axs[0, 2].set_title('(c) 观测面辐照度切片（y2=0）：窗内平顶 ~1 W/m$^2$\n'
                        '竖线 = ±D/2 窗口边缘（原书图 6.8b/6.11b）')
    axs[0, 2].set_xlabel('x2 [mm]'); axs[0, 2].set_ylabel('辐照度 [W/m$^2$]')
    axs[0, 2].legend(fontsize=8); axs[0, 2].grid(alpha=0.3)
    m3 = np.abs(x2mm) <= 8
    an_ph = np.angle(np.exp(1j * (2 * np.pi / wvl / (2 * Dz)
                                  * x2[0, m3] ** 2 - np.pi / 2)))
    axs[1, 0].plot(x2mm[m3], np.angle(U_ob[sl, m3]), 'b.', ms=2,
                   label='numerical')
    axs[1, 0].plot(x2mm[m3], an_ph, 'k-', lw=1.2, label='analytic')
    axs[1, 0].set_title('(d) sinc 模型观测面相位（y2=0）\n'
                        '窗外混叠可见（原书图 6.9）')
    axs[1, 0].set_xlabel('x2 [mm]'); axs[1, 0].set_ylabel('相位 [rad]')
    axs[1, 0].legend(fontsize=8); axs[1, 0].grid(alpha=0.3)
    axs[1, 1].plot(x2mm[m3], np.angle(U_obg[sl, m3]), 'r.', ms=2,
                   label='numerical')
    axs[1, 1].plot(x2mm[m3], an_ph, 'k-', lw=1.2, label='analytic')
    axs[1, 1].set_title('(e) sinc-Gaussian 模型观测面相位（y2=0）\n'
                        '边缘混叠减轻（原书图 6.12）')
    axs[1, 1].set_xlabel('x2 [mm]'); axs[1, 1].set_ylabel('相位 [rad]')
    axs[1, 1].legend(fontsize=8); axs[1, 1].grid(alpha=0.3)
    m4 = np.abs(x2mm) <= 6
    axs[1, 2].plot(x2mm[m4], I_ok[sl, m4] / I_ok.max(), 'b-', lw=1.2,
                   label='修正 (6.92)：质心 (%.2f, %.2f) mm'
                   % (c_ok[0] * 1e3, c_ok[1] * 1e3))
    axs[1, 2].plot(x2mm[m4], I_bad[sl, m4] / I_ok.max(), 'r--', lw=1.2,
                   label='书印 (6.92)：质心 (%.2f, %.2f) mm'
                   % (c_bad[0] * 1e3, c_bad[1] * 1e3))
    axs[1, 2].axvline(rc[0] * 1e3, color='b', ls=':', lw=1)
    axs[1, 2].axvline(-rc[0] * 1e3, color='r', ls=':', lw=1)
    axs[1, 2].set_title('(f) 勘误 (6.92)：离轴窗口 r_c=(1.2, 0.8) mm\n'
                        '书印模型把窗口推到 -r_c（y2=yc 切片）')
    axs[1, 2].set_xlabel('x2 [mm]'); axs[1, 2].set_ylabel('归一化辐照度')
    axs[1, 2].legend(fontsize=8); axs[1, 2].grid(alpha=0.3)
    fig.suptitle('图 6-14（自绘仿真）工程 2：点源模型（sinc / sinc-Gaussian）'
                 '与式 (6.92) 勘误对照（λ=1 μm, Δz=1 m）', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fp = os.path.join(outdir, 'fig-6-14.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
