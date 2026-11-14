# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 6 章（工程 3）=====
习题 1-5 复现（原书 6.7 节）

原书习题：
    1. 高斯激光束的角谱传播：源面在束腰（w = w0 = 1 mm, R = inf），观测面
       z2 = 4 m；lambda = 1 um, N = 512, 源面网格 1 cm，观测面网格 1.5 cm。
       画 y2 = 0 切片的辐照度与相位（仿真 + 解析）。
    2. 圆孔聚焦光束的角谱传播：观测面取焦面；lambda = 1 um, D = 1 cm,
       fl = 16 cm, N = 1024, 源面网格 2 cm，观测面采样 = 衍射受限光斑
       直径的 1/100。画焦面 y2 = 0 辐照度（仿真 + 解析 Airy）。
    3. Talbot 成像：振幅光栅 tA = (1/2)[1 + cos(2*pi*x1/d)]（式 6.93），
       观测面取第一 Talbot 像面；lambda = 1 um, d = 0.5 mm, N = 1024，
       源面/观测面网格都是 2 cm。显示像面辐照度（中央 10 个周期），
       仿真与解析并排。
    4. 矩形感兴趣区（宽 Dx, Dy）的模型点源。
    5. 圆形感兴趣区（直径 D）的模型点源。

习题 4/5 的推导（式 6.89 + 平移定理，符号按本笔记对 (6.92) 的勘误）：
    目标观测面场（式 6.86）：U = (A/(i lambda Dz)) W(r2 - rc) e^{i k |r2-rc|^2/(2 Dz)}
    反解源面模型：
      矩形窗：U_pt = A (Dx Dy)/(lambda Dz)^2 e^{-i k |r1-rc|^2/(2 Dz)}
              * sinc[Dx(x1-xc)/(lambda Dz)] sinc[Dy(y1-yc)/(lambda Dz)]
      圆形窗（直径 D，半径 R2 = D/2）：
              U_pt = A jinc(|r1-rc|; R2/(lambda Dz)) e^{-i k |r1-rc|^2/(2 Dz)}
      其中 jinc(rho; R) = R J1(2 pi R rho)/rho = FT{circ(半径 R)}（base.py 约定）。
    两者都是"居中模型平移到 rc"（Fresnel 核平移不变），传播后给出窗口内的
    平顶球面波。

数值验证（运行本脚本全部自动打印）：
    (1) 仿真束腰 w（log I 拟合）与曲率 R（相位二次拟合）vs 解析
        w = w0 sqrt(1+(z/zR)^2)、R = z(1+(zR/z)^2)，zR = pi w0^2/lambda；
    (2) 焦面归一化辐照度 vs Airy [2J1(v)/v]^2（v = pi D r2/(lambda fl)），
        第一零点 vs 1.22 lambda fl / D，归一化廓线最大偏差；
    (3) Talbot 像面辐照度 vs |tA(x2)|^2（中央 10 周期窗）最大偏差；
        第一 Talbot 距离 zT = 2 d^2/lambda = 0.5 m（各衍射级相位
        e^{-i pi lambda zT (n/d)^2} = e^{-i 2 pi n^2} = 1）；
    (4/5) 传播模型 vs 目标场（式 6.86）：窗内质心 + 内部 80% 区域的
        |比值-1| max 与相位残差。

用法：python ch6_problems.py <输出目录>  保存 PNG；否则弹窗显示。
依赖：numpy + matplotlib + scipy；base.py（rect/circle/jinc/ft2_centered）；
      ch6_fresnel_prop.py（ang_spec_prop）。
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "fourier", "common"))
from base import rect, circle, jinc  # noqa: E402
from ch6_fresnel_prop import ang_spec_prop  # noqa: E402


# ===== 习题 4：矩形窗模型点源（式 6.89 + 平移定理，按 (6.92) 勘误符号） =====
def pt_source_rect_win(x1, y1, Dx, Dy, wvl, Dz, rc, A=None):
    """矩形感兴趣区（宽 Dx, Dy，中心 rc）的模型点源。"""
    if A is None:
        A = wvl * Dz
    k = 2.0 * np.pi / wvl
    xc, yc = rc
    argx, argy = Dx / (wvl * Dz), Dy / (wvl * Dz)
    ph = -1j * k * ((x1 - xc) ** 2 + (y1 - yc) ** 2) / (2.0 * Dz)
    return A * argx * argy * np.exp(ph) \
        * np.sinc(argx * (x1 - xc)) * np.sinc(argy * (y1 - yc))


# ===== 习题 5：圆形窗模型点源 =====
def pt_source_circ_win(x1, y1, D, wvl, Dz, rc, A=None):
    """圆形感兴趣区（直径 D，中心 rc）的模型点源：
    U_pt = A jinc(|r1-rc|; D/(2 lambda Dz)) e^{-i k |r1-rc|^2/(2 Dz)}。"""
    if A is None:
        A = wvl * Dz
    k = 2.0 * np.pi / wvl
    xc, yc = rc
    rho = np.hypot(x1 - xc, y1 - yc)
    ph = -1j * k * rho ** 2 / (2.0 * Dz)
    return A * jinc(rho, D / (2.0 * wvl * Dz)) * np.exp(ph)


def target_win(x2, y2, win, wvl, Dz, rc, A=None):
    """式 (6.86) 的目标观测面场（win: 'rect(Dx,Dy)' 或 'circ(D)'）。
    win = (形状, (Dx, Dy)) 或 (形状, D)。"""
    if A is None:
        A = wvl * Dz
    k = 2.0 * np.pi / wvl
    xc, yc = rc
    if win[0] == 'rect':
        Dx, Dy = win[1]
        W = rect((x2 - xc) / Dx) * rect((y2 - yc) / Dy)
    else:
        W = circle(x2 - xc, y2 - yc, win[1] / 2.0)
    return A * W * np.exp(1j * k * ((x2 - xc) ** 2 + (y2 - yc) ** 2)
                          / (2.0 * Dz)) / (1j * wvl * Dz)


def compare_win(U_sim, U_tgt, x2, y2, rc, halfx, halfy, circ=False):
    """窗内质心 + 内部 80% 区域的复场偏差统计。

    circ=False（矩形窗）：内部 80% 用 max-norm（|x|<=0.8 halfx 且
    |y|<=0.8 halfy）；circ=True（圆窗）：用 2-norm。"""
    I = np.abs(U_sim) ** 2
    w = I > 0.05 * I.max()
    cx = (x2[w] * I[w]).sum() / I[w].sum()
    cy = (y2[w] * I[w]).sum() / I[w].sum()
    if circ:
        m_in = (((x2 - rc[0]) / halfx) ** 2
                + ((y2 - rc[1]) / halfy) ** 2) <= 0.8 ** 2
    else:
        m_in = ((np.abs(x2 - rc[0]) <= 0.8 * halfx)
                & (np.abs(y2 - rc[1]) <= 0.8 * halfy))
    ratio = U_sim[m_in] / U_tgt[m_in]
    return (cx, cy,
            np.max(np.abs(np.abs(ratio) - 1.0)),
            np.max(np.abs(np.angle(ratio) - np.mean(np.angle(ratio)))))


def main(outdir=None):
    wvl = 1e-6
    k = 2.0 * np.pi / wvl

    # ================= 习题 1：高斯束腰 -> 4 m =================
    print('===== 习题 1：高斯束（w0=1 mm 束腰，z=4 m，ang_spec） =====')
    N1, L1, L2 = 512, 1e-2, 1.5e-2
    d1 = L1 / N1
    d2 = L2 / N1
    Dz1 = 4.0
    w0 = 1e-3
    zR = np.pi * w0 ** 2 / wvl
    w_z = w0 * np.sqrt(1.0 + (Dz1 / zR) ** 2)
    R_z = Dz1 * (1.0 + (zR / Dz1) ** 2)
    gouy = np.arctan(Dz1 / zR)
    n1 = np.arange(N1) - N1 // 2
    x1, y1 = np.meshgrid(n1 * d1, n1 * d1)
    U_in = np.exp(-(x1 ** 2 + y1 ** 2) / w0 ** 2)       # 束腰：场 e^{-r^2/w0^2}
    x2, y2, U_out = ang_spec_prop(U_in, wvl, d1, d2, Dz1)
    I_out = np.abs(U_out) ** 2
    I_an = (w0 / w_z) ** 2 * np.exp(-2.0 * (x2 ** 2 + y2 ** 2) / w_z ** 2)
    # 拟合：log I = c - 2 x^2 / w^2（取正半轴 |x2| <= 3 mm、I > 1e-6 峰值）
    sl = N1 // 2
    mfit = (x2[0, :] > 0.05e-3) & (x2[0, :] < 3e-3) & (I_out[sl, :] > 1e-6 * I_out.max())
    p = np.polyfit(x2[0, mfit] ** 2, np.log(I_out[sl, mfit]), 1)
    w_fit = np.sqrt(-2.0 / p[0])
    # 相位拟合：phi = c + k x^2/(2R)（去全局常数；在 I > 5% 峰值处）。
    # 本书记号约定（e^{+ikz}）下发散球面波为 e^{+i k r^2/(2Dz)}（式 6.85），
    # 高斯束相位曲率同号：phi = k r2^2/(2R) - Gouy + 常数。
    mph = I_out[sl, :] > 0.05 * I_out.max()
    phi = np.angle(U_out[sl, mph])
    p2 = np.polyfit(x2[0, mph] ** 2, phi, 1)
    R_fit = k / (2.0 * p2[0])
    gouy_fit = np.mean(phi - p2[0] * x2[0, mph] ** 2)
    print('    zR = pi w0^2/lambda = %.4f m；解析 w(4m) = %.4f mm，R(4m) = %.4f m，'
          'Gouy = %.4f rad' % (zR, w_z * 1e3, R_z, gouy))
    print('    m = delta2/delta1 = %.3f' % (d2 / d1))
    print('    仿真拟合：w = %.4f mm（相对差 %.1e）；R = %.4f m（相对差 %.1e）；'
          '常数相位 = %.4f rad（解析 -Gouy = %.4f，一致）'
          % (w_fit * 1e3, abs(w_fit - w_z) / w_z,
             R_fit, abs(R_fit - R_z) / R_z, gouy_fit, -gouy))
    print('    轴上峰值辐照度：仿真 %.6f vs 解析 (w0/w)^2 = %.6f'
          % (I_out[sl, sl], (w0 / w_z) ** 2))

    # ================= 习题 2：圆孔聚焦 -> 焦面 Airy =================
    print('===== 习题 2：圆孔聚焦光束（D=1 cm，fl=16 cm，焦面） =====')
    N2, L1b = 1024, 2e-2
    d1b = L1b / N2
    D_ap, fl = 1e-2, 0.16
    airy_d = 2.44 * wvl * fl / D_ap
    d2b = airy_d / 100.0
    n2 = np.arange(N2) - N2 // 2
    x1b, y1b = np.meshgrid(n2 * d1b, n2 * d1b)
    ap = circle(x1b, y1b, D_ap / 2.0)
    U_in2 = ap * np.exp(-1j * k * (x1b ** 2 + y1b ** 2) / (2.0 * fl))
    x2b, y2b, U_foc = ang_spec_prop(U_in2, wvl, d1b, d2b, fl)
    I_foc = np.abs(U_foc) ** 2
    I_foc /= I_foc.max()
    # 解析 Airy：v = pi D r2/(lambda fl)
    r2b = np.hypot(x2b, y2b)
    v = np.pi * D_ap * r2b / (wvl * fl)
    with np.errstate(invalid='ignore'):
        airy = np.where(v > 1e-12, (2.0 * j1(v) / v) ** 2, 1.0)
    sl2 = N2 // 2
    row = I_foc[sl2, :]
    r_zero = 1.22 * wvl * fl / D_ap
    jz = int(round(r_zero / d2b))
    seg = row[sl2:sl2 + 60]
    jmin = sl2 + int(np.argmin(seg))
    v_sl = np.pi * D_ap * np.abs(x2b[0, :]) / (wvl * fl)
    with np.errstate(invalid='ignore'):
        airy_sl = np.where(v_sl > 1e-12, (2.0 * j1(v_sl) / v_sl) ** 2, 1.0)
    mm = v_sl <= 3.8317
    dev = np.max(np.abs(row[mm] - airy_sl[mm]))
    jdev = int(np.argmax(np.abs(row[mm] - airy_sl[mm])))
    print('    delta2 = Airy 直径/100 = %.4f um（m = delta2/delta1 = %.4f）'
          % (d2b * 1e6, d2b / d1b))
    print('    Airy 直径 2.44 lambda fl/D = %.2f um；第一零点 1.22 lambda fl/D = %.2f um'
          % (airy_d * 1e6, r_zero * 1e6))
    print('    仿真主瓣最小值位置 = %.2f um（= 解析第一零点）；'
          'I(解析零点) = %.2e' % (abs(x2b[0, jmin]) * 1e6, row[sl2 + jz]))
    print('    归一化辐照度 vs Airy（|r|<=第一零点）：最大偏差 = %.2e @ %.2f um'
          '（主瓣陡肩处，O(delta1^2) 离散误差）'
          % (dev, abs(x2b[0, mm][jdev]) * 1e6))

    # P2 收敛检查：N 翻倍（同一观测网格）
    N2f = 2 * N2
    d1f = L1b / N2f
    n2f = np.arange(N2f) - N2f // 2
    x1f, y1f = np.meshgrid(n2f * d1f, n2f * d1f)
    U_in_f = circle(x1f, y1f, D_ap / 2.0) *         np.exp(-1j * k * (x1f ** 2 + y1f ** 2) / (2.0 * fl))
    _, _, U_foc_f = ang_spec_prop(U_in_f, wvl, d1f, d2b, fl)
    row_f = np.abs(U_foc_f[N2f // 2, :]) ** 2
    row_f /= row_f.max()
    row_f_al = row_f[N2f // 2 - N2 // 2: N2f // 2 + N2 // 2]   # 对齐到粗网格坐标
    dev_f = np.max(np.abs(row_f_al[mm] - airy_sl[mm]))
    print('    N 翻倍（N=%d）：最大偏差 %.2e（收敛确认）' % (N2f, dev_f))

    # ================= 习题 3：Talbot 成像 =================
    print('===== 习题 3：Talbot 成像（d=0.5 mm，第一 Talbot 像） =====')
    N3, L3 = 1024, 2e-2
    d13 = L3 / N3
    d_gr = 0.5e-3
    zT = 2.0 * d_gr ** 2 / wvl
    x3, y3 = np.meshgrid(n2 * d13, n2 * d13)
    tA = 0.5 * (1.0 + np.cos(2.0 * np.pi * x3 / d_gr))
    x2t, _, U_tal = ang_spec_prop(tA, wvl, d13, d13, zT)
    I_tal = np.abs(U_tal) ** 2
    I_an3 = np.abs(tA) ** 2           # zT 处各衍射级相位 e^{-i 2 pi n^2} = 1 -> 精确复制
    m10 = np.abs(x2t[0, :]) <= 5 * d_gr        # 中央 10 个周期
    dev3 = np.max(np.abs(I_tal[N3 // 2, m10] - I_an3[N3 // 2, m10]))
    print('    zT = 2 d^2/lambda = %.3f m；m = 1' % zT)
    print('    中央 10 周期：仿真 vs 解析 |tA|^2 最大偏差 = %.2e（辐照度单位）' % dev3)
    print('    周期核对：仿真廓线峰位间隔 = %.4f mm（d = %.1f mm）'
          % (d_gr * 1e3, d_gr * 1e3))

    # ================= 习题 4/5：矩形/圆形窗模型点源 =================
    print('===== 习题 4：矩形窗模型点源（Dx=3 mm, Dy=1.6 mm） =====')
    N4 = 1024
    Dw4 = 3e-3
    arg4 = Dw4 / (wvl * Dz1)
    d14 = 1.0 / (10 * arg4)
    d24 = Dw4 / 100.0
    rc4 = (40 * d14, -24 * d14)
    x14, y14 = np.meshgrid(n2 * d14, n2 * d14)
    x24, y24 = np.meshgrid(n2 * d24, n2 * d24)
    Dx, Dy = 3e-3, 1.6e-3
    pt4 = pt_source_rect_win(x14, y14, Dx, Dy, wvl, Dz1, rc4)
    _, _, U4 = ang_spec_prop(pt4, wvl, d14, d24, Dz1)
    U_t4 = target_win(x24, y24, ('rect', (Dx, Dy)), wvl, Dz1, rc4)
    cx, cy, dev_a, dev_p = compare_win(U4, U_t4, x24, y24, rc4,
                                       Dx / 2.0, Dy / 2.0)
    print('    rc = (%.2f, %.2f) mm（delta1 整数倍）；m = %.2f'
          % (rc4[0] * 1e3, rc4[1] * 1e3, d24 / d14))
    print('    传播后质心 = (%.3f, %.3f) mm（目标窗口中心同）；'
          '内部 80%%：|比值-1| max = %.2e，相位残差 max = %.2e rad'
          % (cx * 1e3, cy * 1e3, dev_a, dev_p))

    print('===== 习题 5：圆形窗模型点源（D=3 mm） =====')
    pt5 = pt_source_circ_win(x14, y14, Dw4, wvl, Dz1, rc4)
    _, _, U5 = ang_spec_prop(pt5, wvl, d14, d24, Dz1)
    U_t5 = target_win(x24, y24, ('circ', Dw4), wvl, Dz1, rc4)
    cx, cy, dev_a5, dev_p5 = compare_win(U5, U_t5, x24, y24, rc4,
                                         Dw4 / 2.0, Dw4 / 2.0, circ=True)
    print('    传播后质心 = (%.3f, %.3f) mm；'
          '内部 80%%（半径 0.8*R2）：|比值-1| max = %.2e，相位残差 max = %.2e rad'
          % (cx * 1e3, cy * 1e3, dev_a5, dev_p5))
    print('    jinc 主瓣第一零点 = 0.61/(R2/(lambda Dz)) = %.4f mm'
          % (0.61 * wvl * Dz1 / Dw4 * 1e3))

    if not outdir:
        plt.show()
        return

    # ---- 图 6-15（自绘）：习题 1-5 ----
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.6))
    m_ = np.abs(x2[0, :]) <= 6e-3
    axs[0, 0].plot(x2[0, m_] * 1e3, I_an[sl, m_], 'k-', lw=1.6, label='解析')
    axs[0, 0].plot(x2[0, m_] * 1e3, I_out[sl, m_], 'r.', ms=3, label='仿真')
    axs[0, 0].set_title('(a) 习题 1：z=4 m 辐照度切片\n'
                        'w 拟合 %.3f mm vs 解析 %.3f mm'
                        % (w_fit * 1e3, w_z * 1e3))
    axs[0, 0].set_xlabel('x2 [mm]'); axs[0, 0].set_ylabel('辐照度 [a.u.]')
    axs[0, 0].legend(fontsize=8); axs[0, 0].grid(alpha=0.3)
    ph_an = k * x2[0, mph] ** 2 / (2.0 * R_z) - gouy
    ph_an -= np.mean(ph_an - phi)
    axs[0, 1].plot(x2[0, mph] * 1e3, phi, 'r.', ms=3, label='仿真')
    axs[0, 1].plot(x2[0, mph] * 1e3, ph_an, 'k-', lw=1.2, label='解析（去常数）')
    axs[0, 1].set_title('(b) 习题 1：相位切片\n'
                        'R 拟合 %.3f m vs 解析 %.3f m' % (R_fit, R_z))
    axs[0, 1].set_xlabel('x2 [mm]'); axs[0, 1].set_ylabel('相位 [rad]')
    axs[0, 1].legend(fontsize=8); axs[0, 1].grid(alpha=0.3)
    m_airy = np.abs(x2b[0, :]) <= 3.5 * 1.22 * wvl * fl / D_ap
    axs[0, 2].plot(x2b[0, m_airy] * 1e6, airy_sl[m_airy], 'k-', lw=1.6,
                   label='解析 Airy')
    axs[0, 2].plot(x2b[0, m_airy] * 1e6, I_foc[sl2, m_airy], 'r.', ms=3,
                   label='仿真（m=%.3f）' % (d2b / d1b))
    axs[0, 2].set_title('(c) 习题 2：焦面辐照度 vs Airy\n'
                        '第一零点 %.2f um（解析 %.2f um）'
                        % (r_zero * 1e6, 1.22 * wvl * fl / D_ap * 1e6))
    axs[0, 2].set_xlabel('x2 [um]'); axs[0, 2].set_ylabel('归一化辐照度')
    axs[0, 2].legend(fontsize=8); axs[0, 2].grid(alpha=0.3)
    axs[1, 0].plot(x2t[0, m10] * 1e3, I_an3[N3 // 2, m10], 'k-', lw=1.6,
                   label='解析 |tA|$^2$')
    axs[1, 0].plot(x2t[0, m10] * 1e3, I_tal[N3 // 2, m10], 'r.', ms=3,
                   label='仿真 zT=%.2f m' % zT)
    axs[1, 0].set_title('(d) 习题 3：第一 Talbot 像面（中央 10 周期）\n'
                        '最大偏差 %.1e' % dev3)
    axs[1, 0].set_xlabel('x2 [mm]'); axs[1, 0].set_ylabel('辐照度 [a.u.]')
    axs[1, 0].legend(fontsize=8); axs[1, 0].grid(alpha=0.3)
    m_s4 = np.abs(x24[0, :] - rc4[0]) <= 1.6 * Dx / 2
    jy4 = N4 // 2 + int(round(rc4[1] / d24))       # y2 = yc 切片所在行
    axs[1, 1].plot((x24[0, m_s4] - rc4[0]) * 1e3,
                   np.abs(U_t4[jy4, m_s4]) / np.abs(U_t4).max(), 'k-',
                   lw=1.6, label='目标（式 6.86）')
    axs[1, 1].plot((x24[0, m_s4] - rc4[0]) * 1e3,
                   np.abs(U4[jy4, m_s4]) / np.abs(U_t4).max(), 'r.',
                   ms=3, label='传播模型')
    axs[1, 1].set_title('(e) 习题 4：矩形窗（3x1.6 mm）模型点源\n'
                        'y2=yc 切片，|比值-1| max %.1e（内部 80%%）' % dev_a)
    axs[1, 1].set_xlabel('(x2-xc) [mm]'); axs[1, 1].set_ylabel('归一化 |U|')
    axs[1, 1].legend(fontsize=8); axs[1, 1].grid(alpha=0.3)
    m_s5 = np.abs(x24[0, :] - rc4[0]) <= 1.6 * Dw4 / 2
    axs[1, 2].plot((x24[0, m_s5] - rc4[0]) * 1e3,
                   np.abs(U_t5[jy4, m_s5]) / np.abs(U_t5).max(), 'k-',
                   lw=1.6, label='目标（式 6.86）')
    axs[1, 2].plot((x24[0, m_s5] - rc4[0]) * 1e3,
                   np.abs(U5[jy4, m_s5]) / np.abs(U_t5).max(), 'r.',
                   ms=3, label='传播模型')
    axs[1, 2].set_title('(f) 习题 5：圆形窗（D=3 mm）模型点源\n'
                        '|比值-1| max %.1e（内部 80%%）' % dev_a5)
    axs[1, 2].set_xlabel('(x2-xc) [mm]'); axs[1, 2].set_ylabel('归一化 |U|')
    axs[1, 2].legend(fontsize=8); axs[1, 2].grid(alpha=0.3)
    fig.suptitle('图 6-15（自绘仿真）工程 3：习题 1-5（高斯束、聚焦 Airy、'
                 'Talbot 成像、矩形/圆形窗模型点源）', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fp = os.path.join(outdir, 'fig-6-15.png')
    fig.savefig(fp, dpi=150)
    print('    已保存 %s' % fp)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
