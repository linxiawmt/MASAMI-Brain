# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 4 章 =====
透镜的傅里叶变换性质（原书 4.2 节，Listings 4.3-4.5 与图 4.2-4.4）

原书要点：焦距 f_l 的会聚薄镜给场乘相位 exp(-i k/(2 f_l)(x^2+y^2))
（式 4.7）。透明物片放在三个位置、传播到后焦面（Delta_z = f_l）时，
物面上的透镜相位与菲涅耳积分（式 4.6）的二次相位相消，剩下的都是
"常数 prefactor x 二次相位 x 一个 FT"，只是 FT 的频率标度与相位因子不同：
    物贴透镜        式 (4.10)：pref = 1/(i lambda f_l)，相位 k/(2 f_l) r2^2
    物在透镜前 d    式 (4.12)：pref 同上，相位 k/(2 f_l)(1-d/f_l) r2^2，
                    且严格形式带移位光瞳 P(x1 + d/f_l x2)（渐晕，Listing 4.4 未实现）
    物在透镜后 d    式 (4.16)：pref = f_l/d / (i lambda d)，相位 k/(2 d) r2^2，
                    FT 标度 lambda d（几何光学近似，要求 d << f_l）
两个特例：d=0 退回物贴透镜；d=f_l（物在前焦面）时相位因子=1，成为精确 FT。

原书勘误（本脚本按修正后的签名实现）：
    Listing 4.5 第 1-2 行的函数签名 lens_behind_ft(Uin, wvl, d1, f) 少了
    参数 d，而函数体（第 10、16-17 行）用到 d；第 4 行注释里的
    lens_behind_ft(Uin, wvl, d1, d, f) 才是正确签名。另外 Listing 4.4/4.5
    的返回值顺序 [x2 y2 Uout] 与 Listing 4.1/4.3 的 [Uout x2 y2] 不一致
    （能跑，但混用容易翻车）。

数值验证（运行本脚本全部自动打印）：
    (1) d=0 时 lens_in_front_ft 与 lens_against_ft 逐点一致；
    (2) d=f_l 时 lens_in_front_ft 退化为精确 FT（无二次相位），并与
        物贴透镜结果相差恰好 exp(i k/(2 f_l) r2^2)；
    (3) 独立路线对照：直接按菲涅耳衍射积分（式 4.6）对"物面场 =
        透过率 x 透镜相位"做二维数值求和（不用 FFT），与
        lens_against_ft 的 FT 路线逐点对比；
    (4) 能量守恒（离散 Parseval）：sum|U2|^2 d2^2 == 物面积/(lambda f_l)^2；
    (5) 物在透镜后：焦斑宽度 ~ lambda d / D 与峰值 ~ (f_l/d)^2 的实测 vs 理论；
    (6) 渐晕演示：物前 d=0.5 f_l、有限光瞳时，式 (4.11) 的移位光瞳
        P(x1 + d/f_l x2) 使离轴焦点变暗——Listing 4.4（无移位）高估边缘。

MATLAB/Python 差异：
    - 原书 ft2 = fftshift(fft2(fftshift(g)))*delta^2；base.py 的
      ft2_centered 输入端用 ifftshift，偶数 N 下两者逐元素一致。
    - np.meshgrid 默认 'xy'，与 MATLAB meshgrid 一致（x2 沿列变化）。
用法：python ch4_lenses.py <输出目录>  保存 PNG；否则弹窗显示。
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
from base import rect, ft2_centered  # noqa: E402


# ===== 书 Listing 4.3：物贴透镜（lens_against_ft） =====
def lens_against_ft(Uin, wvl, d1, f):
    """物贴透镜 -> 后焦面（式 4.10）。f = f_l 焦距 [m]。"""
    N = Uin.shape[0]
    k = 2 * np.pi / wvl
    fX = (np.arange(N) - N // 2) / (N * d1)
    x2, y2 = np.meshgrid(wvl * f * fX, wvl * f * fX)   # 观察面坐标：标度 lambda*f_l
    Uout = np.exp(1j * k / (2 * f) * (x2 ** 2 + y2 ** 2)) \
        / (1j * wvl * f) * ft2_centered(Uin, d1)
    return Uout, x2, y2


# ===== 书 Listing 4.4：物在透镜前（lens_in_front_ft） =====
def lens_in_front_ft(Uin, wvl, d1, f, d):
    """物在透镜前 d -> 后焦面（式 4.12，取 P=1，即不含移位光瞳的渐晕项）。

    d=0 退化为 lens_against_ft；d=f 时二次相位消失（精确 FT）。
    注意：原书 Listing 4.4 的返回值顺序是 [x2 y2 Uout]，这里统一为
    (Uout, x2, y2)，与其余两个函数保持一致。"""
    N = Uin.shape[0]
    k = 2 * np.pi / wvl
    fX = (np.arange(N) - N // 2) / (N * d1)
    x2, y2 = np.meshgrid(wvl * f * fX, wvl * f * fX)
    Uout = 1 / (1j * wvl * f) \
        * np.exp(1j * k / (2 * f) * (1 - d / f) * (x2 ** 2 + y2 ** 2)) \
        * ft2_centered(Uin, d1)
    return Uout, x2, y2


# ===== 书 Listing 4.5：物在透镜后（lens_behind_ft，勘误修正签名） =====
def lens_behind_ft(Uin, wvl, d1, d, f):
    """物在透镜后（距焦面 d）-> 焦面（式 4.16）。

    原书 Listing 4.5 的签名缺 d（勘误），这里补上；FT 标度是 lambda*d，
    幅度乘 f_l/d（会聚球面波在物面处的振幅，式 4.13）。"""
    N = Uin.shape[0]
    k = 2 * np.pi / wvl
    fX = (np.arange(N) - N // 2) / (N * d1)
    x2, y2 = np.meshgrid(wvl * d * fX, wvl * d * fX)   # 标度是 lambda*d ！
    Uout = f / d * 1 / (1j * wvl * d) \
        * np.exp(1j * k / (2 * d) * (x2 ** 2 + y2 ** 2)) \
        * ft2_centered(Uin, d1)
    return Uout, x2, y2


def fresnel_direct_slice(U1, wvl, d1, Dz, x2_vals, y2=0.0):
    """独立路线：按菲涅耳衍射积分（式 4.6）对全部网格点直接二维求和
    （不用 FFT），返回 y=y2 行各 x2_vals 处的 U。

    U(x2,y2) = e^{ik Dz}/(i lambda Dz) sum U1(x1,y1)
               x e^{i k/2Dz (x1^2+y1^2)} e^{-i 2pi/(lambda Dz)(x2 x1 + y2 y1)} d1^2
    （e^{ik Dz} 轴上相位与原书一样略去；求和点数 = N^2，演示用小网格。）"""
    N = U1.shape[0]
    k = 2 * np.pi / wvl
    n = np.arange(N) - N // 2
    xa = n * d1
    X1, Y1 = np.meshgrid(xa, xa)
    quad = np.exp(1j * k / (2 * Dz) * (X1 ** 2 + Y1 ** 2))    # 源面二次相位
    kern = np.exp(-2j * np.pi / (wvl * Dz)
                  * (x2_vals[:, None] * X1.ravel()[None, :]
                     + y2 * Y1.ravel()[None, :]))
    src = (U1 * quad).ravel()
    return np.einsum('j,ij->i', src, kern) * d1 ** 2 / (1j * wvl * Dz)


def main(outdir=None):
    # ---- 演示参数：0.5 mm 方孔物片，f_l = 1 m，lambda = 1 um ----
    N = 256
    L = 2.5e-3         # 源面网格宽度 [m]
    d1 = L / N
    D = 0.5e-3         # 方孔边长 [m]
    wvl = 1e-6
    k = 2 * np.pi / wvl
    fl = 1.0           # 焦距 [m]
    n = np.arange(N) - N // 2
    x1, y1 = np.meshgrid(n * d1, n * d1)
    tA = rect(x1 / D) * rect(y1 / D)
    area = D ** 2

    print('===== 验证 1：d=0 时 lens_in_front_ft == lens_against_ft =====')
    U_ag, x2, y2 = lens_against_ft(tA, wvl, d1, fl)
    U_bf0, _, _ = lens_in_front_ft(tA, wvl, d1, fl, 0.0)
    print('    max|U_ag - U_bf(d=0)| = %.2e（相对 %.2e）'
          % (np.max(np.abs(U_ag - U_bf0)),
             np.max(np.abs(U_ag - U_bf0)) / np.max(np.abs(U_ag))))

    print('===== 验证 2：d=f_l 时退化为精确 FT；与贴透镜差一个二次相位 =====')
    U_bf1, _, _ = lens_in_front_ft(tA, wvl, d1, fl, fl)
    U_exact = 1 / (1j * wvl * fl) * ft2_centered(tA, d1)   # 纯 FT，无相位因子
    ph_dev = np.abs(np.angle(U_bf1 / U_exact))
    print('    max|angle(U_bf(d=f_l)/U_exact)| = %.2e rad（精确 FT，零相位）'
          % ph_dev.max())
    U_exact_ag = U_exact * np.exp(1j * k / (2 * fl) * (x2 ** 2 + y2 ** 2))
    print('    max|U_ag - exp(ik/2f_l r2^2) * U_exact| = %.2e'
          % np.max(np.abs(U_ag - U_exact_ag)))
    print('    焦面二次相位不小儿科：贴透镜时 r2=1 cm 处相位 = %.1f rad'
          % (k / (2 * fl) * 0.01 ** 2))

    print('===== 验证 3：两条路线对照（直接积分 vs FT，物贴透镜） =====')
    # 物面场 = tA * exp(-i k/2f_l (x1^2+y1^2))（透镜紧贴物片之后），
    # 用式 (4.6) 直接传播 Delta_z = f_l；透镜相位与积分的二次相位相消。
    U1 = tA * np.exp(-1j * k / (2 * fl) * (x1 ** 2 + y1 ** 2))
    sl = N // 2
    xs = x2[sl, :]
    U_dir = fresnel_direct_slice(U1, wvl, d1, fl, xs)
    rel = np.max(np.abs(np.abs(U_dir) - np.abs(U_ag[sl, :]))) \
        / np.max(np.abs(U_ag[sl, :]))
    print('    y2=0 切片：max||U_direct|-|U_ft||/max|U_ft| = %.2e' % rel)

    print('===== 验证 4：能量守恒（离散 Parseval） =====')
    d2 = wvl * fl / L
    P_num = np.sum(np.abs(U_ag) ** 2) * d2 ** 2
    print('    sum|U2|^2 d2^2 = %.5e，物面积 = %.5e（相对偏差 %.2e，'
          '来自节点采样把孔径边界样本整格计入/漏计）'
          % (P_num, area, abs(P_num - area) / area))
    print('    轴上辐照度（物贴透镜）= %.6f W/m^2，理论 (D^2/(lambda f_l))^2 = %.6f'
          % (np.abs(U_ag[sl, sl]) ** 2, (area / (wvl * fl)) ** 2))

    print('===== 验证 5：物在透镜后——焦斑宽度与峰值随 d 变化 =====')
    for d in [0.5 * fl, 0.25 * fl]:
        U_bh, x2b, y2b = lens_behind_ft(tA, wvl, d1, d, fl)
        Ib = np.abs(U_bh) ** 2
        row = Ib[N // 2, :]
        peak_i = row.argmax()
        # 实测主瓣半径：峰值右侧第一个 <1% 峰值的样本（sinc^2 主瓣内
        # 只有 >=1% 峰值的台座，sinc(u)^2=0.01 在 u=0.988，误差 <1.2%）
        right = peak_i + np.argmax(row[peak_i:] < 0.01 * row[peak_i])
        x2_axis = x2b[N // 2, :]
        mainlobe = 2 * (x2_axis[right] - x2_axis[peak_i])
        lobe_th = 2 * wvl * d / D          # sinc^2 主瓣全宽 = 2 lambda d / D
        peak_th = (fl / d) ** 2 * (area / (wvl * d)) ** 2
        print('    d=%.2f m：实测峰值 %.4e vs 理论 %.4e；'
              '实测主瓣全宽 %.3e m vs 理论 %.3e m'
              % (d, row.max(), peak_th, mainlobe, lobe_th))

    print('===== 验证 6：渐晕——式 (4.11) 的移位光瞳 vs Listing 4.4 =====')
    # 物前 d=0.5 f_l，有限圆瞳 D_p=2 mm。Listing 4.4 算的是无移位光瞳
    # （相当于 P=1 或"光瞳跟物片重合"）；严格式 (4.11) 对每个 x2 用
    # P(x1 + d/f_l x2)。这里对 y2=0 切片逐点直接求和比较两者。
    d = 0.5 * fl
    Dp = 2e-3
    P0 = circle_p(x1, y1, Dp)                      # 物面处的光瞳（无移位）
    U_no = lens_in_front_ft(tA * P0, wvl, d1, fl, d)
    # 严格版：直接求和 sum tA(x1) P(x1 + d/f_l x2) e^{-i2pi/(lambda f_l) x2 x1}
    inside = tA != 0
    xa, ya = x1[inside], y1[inside]
    xs6 = xs[np.abs(xs) <= 2.5e-3]
    U_strict = np.empty(xs6.size, dtype=complex)
    for j, xv in enumerate(xs6):
        sh = d / fl * xv
        Pv = circle_p(xa + sh, ya, Dp)             # 移位光瞳
        U_strict[j] = np.sum(tA[inside] * Pv
                             * np.exp(-2j * np.pi / (wvl * fl) * xv * xa)) \
            * d1 ** 2 / (1j * wvl * fl) \
            * np.exp(1j * k / (2 * fl) * (1 - d / fl) * xv ** 2)
    U_no_sl = np.abs(U_no[0][sl, :])[np.abs(xs) <= 2.5e-3]
    for xv in [0.0, 1.0e-3, 1.5e-3, 2.0e-3]:
        j = np.argmin(np.abs(xs6 - xv))
        ratio = np.abs(U_strict[j]) / np.abs(U_no_sl[j])
        print('    x2 = %.1f mm：|U_严格(移位光瞳)|/|U_Listing4.4| = %.3f'
              % (xv * 1e3, ratio))
    print('    （主瓣内 |x2|<=1 mm 两者一致——渐晕尚未发生；主瓣外移位光瞳'
          '截断填掉了 sinc 零点，图样定性不同。Listing 4.4 只在 P=1 或'
          '光瞳随物片一起移动时严格成立）')

    if outdir:
        # ---- 图 fig-4-7：三种物位置的焦面辐照度与相位 ----
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
        U_bf05, _, _ = lens_in_front_ft(tA, wvl, d1, fl, 0.5 * fl)
        U_bh, x2b, _ = lens_behind_ft(tA, wvl, d1, 0.5 * fl, fl)
        m5 = np.abs(x2[sl, :]) <= 2.5e-3
        axes[0].plot(x2[sl, :][m5] * 1e3, np.abs(U_ag[sl, :][m5]) ** 2,
                     'k-', lw=1.6, label='物贴透镜（=物前 d=0）')
        axes[0].plot(x2[sl, :][m5] * 1e3, np.abs(U_bf05[sl, :][m5]) ** 2,
                     'b.', ms=5, label='物前 d=0.5$f_l$（辐照度重合）')
        m6 = np.abs(x2b[sl, :]) <= 2.5e-3
        axes[0].plot(x2b[sl, :][m6] * 1e3, np.abs(U_bh[sl, :][m6]) ** 2,
                     'r--', lw=1.4, label='物后 d=0.5$f_l$（标度 $\\lambda d$）')
        axes[0].set_xlabel(r'$x_2$ [mm]')
        axes[0].set_ylabel(r'辐照度 [W/m$^2$]')
        axes[0].set_title('(a) 焦面辐照度：物前位置不影响图样，物后改变标度')
        axes[0].legend(fontsize=9)
        m7 = np.abs(x2[sl, :]) <= 2e-3
        x7 = x2[sl, :][m7]
        x_fine = np.linspace(-2e-3, 2e-3, 2001)
        for d, col, lab in [(0.0, 'k', '物贴透镜（曲率 $1/f_l$）'),
                            (0.5 * fl, 'tab:blue', '物前 d=0.5$f_l$（$0.5/f_l$）'),
                            (fl, 'tab:green', '物前 $d=f_l$（精确 FT，零相位）')]:
            phi = np.angle(np.exp(1j * k / (2 * fl) * (1 - d / fl)
                                  * x_fine ** 2))
            axes[1].plot(x_fine * 1e3, phi, color=col, lw=1.2, label=lab)
        axes[1].plot(x7 * 1e3, np.angle(U_ag[sl, :][m7]), 'k.', ms=6)
        axes[1].plot(x7 * 1e3, np.angle(U_bf05[sl, :][m7]), 'b.', ms=6)
        axes[1].set_xlabel(r'$x_2$ [mm]')
        axes[1].set_ylabel('相位 [rad]')
        axes[1].set_title('(b) 焦面相位记得物的位置（点=数值，线=理论包裹相位）')
        axes[1].legend(fontsize=8, loc='upper center')
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.tight_layout()
        fp = os.path.join(outdir, 'fig-4-7.png')
        fig.savefig(fp, dpi=150)
        print('\n已保存：%s' % fp)
    else:
        plt.show()


def circle_p(X, Y, D):
    """直径 D 的圆瞳（|r| <= D/2 取 1）。"""
    return np.where(np.hypot(X, Y) <= D / 2, 1.0, 0.0)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
