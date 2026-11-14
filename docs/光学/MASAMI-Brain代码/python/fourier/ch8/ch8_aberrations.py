# -*- coding: utf-8 -*-
"""第 8 章波前像差函数库:seidel_5 / 像差瞳函数 / PSF-MTF 工具 / Zernike(书 8.2-8.3 节)。

与 matlab/fourier/ch8/seidel_5.m 逐行对应;瞳函数按书式 (8.3) 实现,
PSF/MTF 按书式 (8.7)/(8.9) 与 lens_psfmtf 脚本的 FFT 口径实现,供
ch8_lens_psfmtf、ch8_psf_map、ch8_image_super、ch8_exercises 导入。

import 路径(本章脚本统一写法):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "common"))
    from base import coords, circle, rect          # 公共积木(ch1 起)
    from ch8_aberrations import seidel_5, psf_from_H, mtf_from_psf

MATLAB/Python 差异(详见笔记第 8 章坑与陷阱):
  - numpy 索引 0-based(MATLAB 1-based):居中排列的"零坐标"样本,
    MATLAB 是 (M/2+1, M/2+1),numpy 是 (M//2, M//2);剖面取 h2[M//2, :]。
  - MATLAB circshift(a, [r, c]) 的第一变量是行(= y/v 方向),
    numpy 对应 np.roll(a, (r, c), axis=(0, 1)),轴顺序相同、索引基不同。
  - MATLAB nthroot(I, n) 对非负数即 I**(1/n);负底数时 Python ** 返回复数,
    本章辐照度均非负,可直接用 **。
  - FFT 排布:脚本先把频率坐标 fu=fftshift(fu) 预移位,循环里 H 天然处于
    "FFT 序"(DC 在 [0,0]),ifft2 后只需 ifftshift 回中(lens_psfmtf 则是
    居中坐标 + fftshift(H),两种写法等价,书 8.7 节特别说明前者是省时技巧)。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
from base import circle, rect                   # noqa: E402  公共积木(ch1 起)


# ===== 书 8.2.2 节:seidel_5(与 matlab/fourier/ch8/seidel_5.m 逐行对应) =====
def seidel_5(u0, v0, X, Y, wd, w040, w131, w222, w220, w311):
    """前 5 项 Seidel 像差 + 离焦的波前 OPD(书式 (8.2))。

    u0, v0:归一化像面坐标(标量);X, Y:归一化出瞳坐标网格(meshgrid)。
    系数单位为长度(米);书习惯用"波长数"乘 lambda 传入。
    先把瞳坐标旋转 beta = atan2(v0, u0),使 x' 轴指向像点方向
    (旋转对称系统只需在子午面内写多项式,书 8.2.1 节)。
    """
    beta = np.arctan2(v0, u0)                           # 书 15 行:image rotation angle
    u0r = np.sqrt(u0 ** 2 + v0 ** 2)                    # 书 16 行:image height
    # 书 19-20 行:rotate grid
    Xr = X * np.cos(beta) + Y * np.sin(beta)
    Yr = -X * np.sin(beta) + Y * np.cos(beta)
    # 书 23-29 行:Seidel polynomials(式 8.2)
    rho2 = Xr ** 2 + Yr ** 2
    w = (wd * rho2
         + w040 * rho2 ** 2
         + w131 * u0r * rho2 * Xr
         + w222 * u0r ** 2 * Xr ** 2
         + w220 * u0r ** 2 * rho2
         + w311 * u0r ** 3 * Xr)
    return w


# ===== 书 8.3.1 节:像差瞳函数(式 8.3) =====
def pupil_aberrated(fu_grid, fv_grid, lam, zxp, wxp, u0, v0, coeffs):
    """出瞳处的相干传递函数 H(式 8.3 + 8.4 合并)。

    fu_grid, fv_grid:像面频率坐标(cyc/m,任意排布——逐点计算,与排布无关);
    物理瞳坐标 x_p = -lambda*zxp*f_U,归一化 xhat = x_p/wxp = -2*lambda*fnum*f_U
    (两种写法相同,因为 wxp = zxp/(2*fnum));coeffs = (wd,w040,w131,w222,w220,w311)。
    返回 H = circ(f_r/f0) * exp(-j k W),f0 = wxp/(lam*zxp) 为相干截止频率。
    """
    k = 2 * np.pi / lam
    xhat = -lam * zxp * fu_grid / wxp                   # 式 (8.4) 宗量(归一化瞳坐标)
    yhat = -lam * zxp * fv_grid / wxp
    W = seidel_5(u0, v0, xhat, yhat, *coeffs)
    f0 = wxp / (lam * zxp)
    # circ(f_r / f0):Goodman 约定(分母=半径),公共模块 base.circle(X,Y,R)
    # 即"半径 R 圆盘",故传 f0(ch1 勘误:与表 1.2 印刷定义差 2 倍半径)
    return circle(fu_grid, fv_grid, f0) * np.exp(-1j * k * W)


def psf_from_H(H, pre_shifted):
    """由相干传递函数得 PSF(式 8.7):h2 = |IFFT{H}|^2。

    pre_shifted=True:H 已处于 FFT 序(fu 曾做 fftshift 预移位,psf_map 口径),
      直接 ifft2 再 ifftshift 回中;
    pre_shifted=False:H 处于居中排列(lens_psfmtf 口径),
      需 ifftshift(H) -> ifft2 -> fftshift 回中(numpy 与 MATLAB 的
      fftshift/ifftshift 在偶数 N 下互为同一操作)。
    """
    if pre_shifted:
        h2 = np.abs(np.fft.ifftshift(np.fft.ifft2(H))) ** 2
    else:
        h2 = np.fft.fftshift(np.abs(np.fft.ifft2(np.fft.ifftshift(H))) ** 2)
    return h2


def mtf_from_psf(h2):
    """由 PSF 得 MTF(式 8.9,lens_psfmtf 54-56 行口径):|FFT{h2}| 归一化 DC=1。

    h2 为居中排列;内部 ifftshift -> fft2 -> 取模 -> 除 DC -> ifftshift 回中。
    取模之后 fftshift 与否只影响排布,不影响数值。"""
    MTF = np.abs(np.fft.fft2(np.fft.ifftshift(h2)))
    MTF = MTF / MTF[0, 0]                               # 书 55 行:normalize DC to 1
    return np.fft.ifftshift(MTF)


def diffraction_limited_mtf_1d(fu, twof0):
    """一维衍射限 MTF(书式 7.31,lens_psfmtf 59-61 行):圆瞳自相关解析式。"""
    nu = fu / twof0
    mtf = (2 / np.pi) * (np.arccos(np.clip(nu, -1, 1))
                         - nu * np.sqrt(np.clip(1 - nu ** 2, 0, None)))
    return mtf * rect(fu / (2 * twof0))                 # 书 61 行:zero after cutoff


# ===== 习题 8.6:Zernike 标准多项式(第 1-11 项,表 8.4 的多项式列) =====
def zernike_11(X, Y, coeffs):
    """按表 8.4 的 ZEMAX"standard"Zernike 形式逐项列出(与 seidel_5 同风格)。

    p = 归一化极径,A = atan2(Y, X)(相对 x 轴,书 8.8 节习题说明);
    coeffs:长度 >=11 的系数序列(单位:米;表 8.4 按"波长数"给出乘 lambda)。
    第 1 项是 piston(piston 不影响 PSF/MTF);像点位置不显式出现——
    这是 Zernike 与 Seidel 的关键区别(书习题 8.6)。
    """
    c = list(coeffs) + [0.0] * (11 - len(coeffs))
    p = np.hypot(X, Y)
    A = np.arctan2(Y, X)
    Z = [
        c[0] * np.ones_like(X),                                   # 1  piston
        c[1] * np.sqrt(4) * p * np.cos(A),                        # 2  tilt (cos)
        c[2] * np.sqrt(4) * p * np.sin(A),                        # 3  tilt (sin)
        c[3] * np.sqrt(3) * (2 * p ** 2 - 1),                     # 4  defocus
        c[4] * np.sqrt(6) * p ** 2 * np.sin(2 * A),               # 5  astig 45
        c[5] * np.sqrt(6) * p ** 2 * np.cos(2 * A),               # 6  astig 0
        c[6] * np.sqrt(8) * (3 * p ** 3 - 2 * p) * np.sin(A),     # 7  coma (sin)
        c[7] * np.sqrt(8) * (3 * p ** 3 - 2 * p) * np.cos(A),     # 8  coma (cos)
        c[8] * np.sqrt(8) * p ** 3 * np.sin(3 * A),               # 9  trefoil (sin)
        c[9] * np.sqrt(8) * p ** 3 * np.cos(3 * A),               # 10 trefoil (cos)
        c[10] * np.sqrt(5) * (6 * p ** 4 - 6 * p ** 2 + 1),       # 11 spherical
    ]
    w = np.zeros_like(X)
    for z in Z:
        w = w + z
    return w


# ===== 书 8.7.2 节:合成 USAF-1951 风格测试靶(替代书代码分发的 USAF1951B250.png) =====
def make_usaf_target(M, fbase=16.0):
    """合成 USAF 1951 风格分辨力靶(M x M,取值 0..1)。

    原书 image_super 读入随书代码分发的 USAF1951B250.png;本笔记 Python 侧
    用合成靶替代:两列三条靶,元素 n 的频率 = fbase*2^((6-n)/6)(周期/幅宽),
    条宽 = 半周期(占空比 50%),条长 = 4 倍条宽(细长比模仿真实靶)。
    布局:左列竖条、右列横条、中央嵌套方块。坐标归一化到 [0,1],按
    "行 0 = 底"(imshow origin='lower')直接画出,无需书第 4 行的 flipud。
    """
    img = np.zeros((M, M))
    jj, ii = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    xn, yn = ii / M, jj / M                             # xn:列方向(u),yn:行方向(v)

    def tri_bars(x0, y0, freq, vertical, L=0.10):
        """在 (x0,y0) 画三条;vertical=True 条沿 v 方向长、沿 u 排开。"""
        period = 1.0 / freq
        w = period / 2
        for b in range(3):
            if vertical:
                sel = ((xn >= x0 + b * period) & (xn < x0 + b * period + w)
                       & (yn >= y0) & (yn < y0 + L))
            else:
                sel = ((xn >= x0) & (xn < x0 + L)
                       & (yn >= y0 + b * period) & (yn < y0 + b * period + w))
            img[sel] = 1.0

    # 左列:竖直三条组(自下而上频率递减 -> 条越来越大;条长 = 4 倍条宽)
    for n in range(1, 7):
        freq = fbase * 2 ** ((6 - n) / 6)
        period = 1.0 / freq
        y0 = 0.05 + (n - 1) * (4 * period + 0.02)
        tri_bars(0.07, y0, freq, vertical=True, L=4 * period)
    # 右列:水平三条组
    for n in range(1, 7):
        freq = fbase * 2 ** ((6 - n) / 6)
        period = 1.0 / freq
        x0 = 0.60 + 0.03
        y0 = 0.05 + (n - 1) * (4 * period + 0.02)
        tri_bars(x0, y0, freq, vertical=False, L=4 * period)
    # 中央嵌套方块(粗结构,模仿书图中央)
    box(img, 0.40, 0.42, 0.14, 0.14, 1.0)
    box(img, 0.435, 0.455, 0.07, 0.07, 0.0)
    box(img, 0.4525, 0.4725, 0.035, 0.035, 1.0)
    return img


def box(img, x0, y0, w, h, val):
    """实心方块辅助(归一化坐标,左下原点)。"""
    M = img.shape[0]
    jj, ii = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    xn, yn = ii / M, jj / M
    img[(xn >= x0) & (xn < x0 + w) & (yn >= y0) & (yn < y0 + h)] = val
