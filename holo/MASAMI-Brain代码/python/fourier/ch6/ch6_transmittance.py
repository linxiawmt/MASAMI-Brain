# -*- coding: utf-8 -*-
"""第 6 章透过率函数库:tilt / focus / 透镜 / 光栅周期结构(书 6.1-6.4 节)。

与 matlab/fourier/ch6/tilt.m、focus.m 逐行对应;透镜与光栅的透过率按
书式 (6.12)/(6.20)/(6.25) 实现,供本章各脚本(ch6_sqr_beam_*、
ch6_lens_airy、ch6_gratings、ch6_exercises)导入。

import 路径(本章脚本统一写法):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "common"))
    from base import coords, circle
    from ch6_transmittance import tilt, focus

MATLAB/Python 差异(详见笔记第 6 章坑与陷阱):
  - numpy 索引 0-based(MATLAB 1-based):居中排列的"零坐标"样本,
    MATLAB 是 (M/2+1, M/2+1),numpy 是 (M//2, M//2)。
  - 书式 (6.13)/(6.29) 的 circ 按附录 C.1 / Goodman 约定(半径 = 分母),
    与表 1.2 的印刷定义差 2 倍半径(全书约定冲突,ch1 已标勘误)。
    本模块统一用公共模块 base.circle(X, Y, R)(半径 R 的圆盘)实现:
    透镜孔径半径 wL 写 circle(X, Y, wL),波带片半径 w 写 circle(X, Y, w)。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
from base import circle                        # noqa: E402  公共积木(ch1 起)


# ===== 书 6.1 节:tilt(与 matlab/fourier/ch6/tilt.m 逐行对应) =====
def tilt(uin, L, lam, alpha, theta):
    """倾斜透过率(式 6.3):t = exp[jk(x cosθ + y sinθ) tanα]。

    给波前加线性相位斜面,光束朝 (cosθ, sinθ) 方向偏转 α,
    瞄准点 r = z tanα(图 6.2)。倾斜把源频谱平移 tanα/λ(式 6.4/6.5),
    因此受带宽判据 α ≤ λ(1/(2Δx1) − B1)(式 6.6)限制。
    """
    M, N = uin.shape                                    # 书 11 行
    dx = L / M                                          # 书 12 行
    k = 2 * np.pi / lam                                 # 书 13 行
    x, _ = coords_local(L, M)                           # 书 15 行:空间坐标(居中)
    X, Y = np.meshgrid(x, x)                            # 书 16 行
    uout = uin * np.exp(1j * k * (X * np.cos(theta) + Y * np.sin(theta))
                        * np.tan(alpha))                # 书 18-19 行:apply tilt
    return uout


# ===== 书 6.2 节:focus(与 matlab/fourier/ch6/focus.m 逐行对应) =====
def focus(uin, L, lam, zf):
    """聚焦/发散透过率(式 6.9):t = exp[-j(k/2zf)(x²+y²)]。

    zf > 0 会聚(实焦点),zf < 0 发散(虚焦点)。
    与菲涅耳脉冲响应 h(式 4.23)同形的 chirp,指数符号为负。
    聚焦加宽源带宽约 D1/(2λ|zf|)(式 6.10),故有式 (6.11) 的 |zf| 下界。
    """
    M, N = uin.shape                                    # 书 10 行
    dx = L / M                                          # 书 11 行
    k = 2 * np.pi / lam                                 # 书 12 行
    x, _ = coords_local(L, M)                           # 书 14 行
    X, Y = np.meshgrid(x, x)                            # 书 15 行
    uout = uin * np.exp(-1j * k / (2 * zf) * (X ** 2 + Y ** 2))   # 书 17 行
    return uout


# ===== 习题 6.3(a):柱面聚焦(chirp 只作用在 x 轴) =====
def focus_cyl(uin, L, lam, zf):
    """柱面透镜聚焦:exp[-j(k/2zf)x²],y 方向不变(习题 6.3)。"""
    M, N = uin.shape
    k = 2 * np.pi / lam
    x, _ = coords_local(L, M)
    X, _ = np.meshgrid(x, x)     # meshgrid 第一输出:沿列(第 2 维)变化的 x 坐标
    #            (u 的约定:行 = y、列 = x;MATLAB 侧用 [X,~]=meshgrid(x,x))
    return uin * np.exp(-1j * k / (2 * zf) * X ** 2)


# ===== 书 6.3 节:理想薄透镜透过率(式 6.12/6.13) =====
def lens_ta(X, Y, lam, f, wL=None):
    """t_A(x,y) = P(x,y)·exp[-j(k/2f)(x²+y²)](式 6.12)。

    P 为瞳函数(式 6.13):半径 wL 的圆孔。注意式 (6.13) 的 circ 宗量
    是 √(x²+y²)/wL(Goodman 约定,半径 = 分母 = wL);base.circle(X,Y,R)
    画的是半径 R 圆盘,故传 R = wL。
    wL 为 None 时不加瞳(纯相位,无限大透镜)。
    """
    k = 2 * np.pi / lam
    t = np.exp(-1j * k / (2 * f) * (X ** 2 + Y ** 2))
    if wL is not None:
        t = t * circle(X, Y, wL)          # 瞳函数:circ(√(x²+y²)/wL),半径 wL
    return t


# ===== 习题 6.5:波带片(式 6.29 / 6.30) =====
def zone_plate(X, Y, lam, f, w, binary=False):
    """波带片透过率,t_A(r) = ½[1+cos(πr²/(λf))]·circ(r/w)(式 6.29)。

    注:πr²/(λf) = (k/2f)r²。binary=True 时换二值型(式 6.30):
    ½[1+sign(cos(πr²/(λf)))]·circ(r/w)——透过率只有 1/0,可用
    悬浮不透明材料环带做成(反射光栅的电脑模型同理,习题 6.7)。
    circ 约定同 lens_ta:半径 = w(附录 C.1 的 Goodman 约定)。
    """
    k = 2 * np.pi / lam
    r2 = X ** 2 + Y ** 2
    c = np.cos(0.5 * k / f * r2)                       # cos(πr²/(λf))
    if binary:
        t = 0.5 * (1.0 + np.sign(c))
    else:
        t = 0.5 * (1.0 + c)
    return t * circle(X, Y, w)                         # circ(r/w),半径 w


# ===== 附录 C.5:ucomb(与 matlab/fourier/ch6/ucomb.m 逐行对应) =====
def ucomb(x):
    """单位样本"梳":x 为整数处取 1,先四舍五入到 10^-6 截断舍入误差。"""
    x = np.round(np.asarray(x, dtype=float) * 1e6) / 1e6   # 书:round to 10^6ths
    return np.where(np.remainder(x, 1) == 0, 1.0, 0.0)     # 书:rem(x,1)==0


# ===== 书 6.4.2 节:方波光栅的一维构造(与书脚本第 15-19 行对应) =====
def square_wave_1d(x1, P):
    """rect(x/(P/2)) ⊛ (1/P)comb(x/P) 的一维方波(式 6.25 的卷积部分)。

    按书做法用 FFT 求卷积:fc=FFT{ucomb(x/P)},fr=FFT{rect(x/(P/2))},
    ux=IFFT{fc.*fr};不需要额外定标——(a) ifft 已为一次 fft 正确
    定标;(b) 每个单位样本 delta 像有 1/dx1 的面积,抵消第二个 fft 的 dx1。
    前提:P/dx1 必须是整数(ucomb 才能落位),P/(2dx1) 为奇数与 rect 定义一致。
    """
    fc = np.fft.fft(np.fft.fftshift(ucomb(x1 / P)))
    fr = np.fft.fft(np.fft.fftshift(rect_1d(x1 / (P / 2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc * fr))            # 1D conv rect & comb
    return np.real(ux)          # 方波是实函数,舍掉 1e-16 级虚部便于观察


# ===== 习题 6.7(c)/6.8:三角波光栅的一维构造 =====
def tri_wave_1d(x1, P):
    """Λ(x/(P/2)) ⊛ (1/P)comb(x/P) 的一维三角波(式 6.32 的卷积部分)。

    F{Λ(x/(P/2))} = (P/2)sinc²(Pf/2):阶的权重 ∝ sinc²(n/2) ≥ 0,
    与方波的 sinc(n/2)(可正可负)不同——这就是"闪耀"式的相位结构。"""
    fc = np.fft.fft(np.fft.fftshift(ucomb(x1 / P)))
    ft = np.fft.fft(np.fft.fftshift(tri_1d(x1 / (P / 2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc * ft))
    return np.real(ux)


# ---------- 以下为附录 C 的 rect/tri(numpy 版,符号处理) ----------
def rect_1d(x):
    """rect(x):|x|<=1/2 取 1(附录 C.3 的 out=abs(x)<=1/2)。"""
    return np.where(np.abs(x) <= 0.5, 1.0, 0.0)


def tri_1d(x):
    """Λ(x):|x|<=1 取 1-|x|(附录 C.4)。"""
    return np.where(np.abs(x) <= 1.0, 1.0 - np.abs(x), 0.0)


def coords_local(L, M):
    """书 2.4 节坐标(居中):x = -L/2 : dx : L/2-dx。

    直接复用公共模块 base.coords 的空间轴(第 5 章起的约定)。
    """
    from base import coords
    return coords(L, M)
