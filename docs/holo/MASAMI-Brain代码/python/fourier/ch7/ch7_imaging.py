# -*- coding: utf-8 -*-
"""第 7 章成像仿真核心模块(书 7.2/7.3 节)。

把书 coh_image / incoh_image 两个脚本里"可复用"的部分抽成函数:
  - 相干传递函数 H(式 7.7):方瞳(式 7.9)、圆瞳(式 7.13)、矩形瞳、环瞳
  - 相干成像管线(书 coh_image 第 35-39 行,式 7.16)
  - OTF 的 FFT 自相关算法(书 incoh_image 第 22-23 行,式 7.26)
  - 非相干成像管线(书 incoh_image 第 32-36 行,式 7.33)
  - 圆瞳 OTF 解析式(式 7.31)与方瞳 OTF 解析式(式 7.29)
  - 附录 C.6 的 udelta(习题 7.10 用)

MATLAB/Python 差异(详见笔记第 7 章坑与陷阱):
  - 阵列"排列"约定:u1/Ig/H 按居中排列存放;进 FFT 前 fftshift、
    出来后 ifftshift(偶数 N 时两者相同)。书 incoh_image 的 OTF 故意
    留在 FFT 排列(第 22-23 行不 fftshift),与 fft2(fftshift(Ig)) 对乘
    ——这里原样保留,注释里标明。
  - MATLAB 的 OTF(1,1) 是 1-based 的 (1,1) = FFT 排列的零频样本;
    numpy 0-based 对应 OTF[0, 0]。
  - base.circle(X, Y, R) 画半径 R 的圆盘,等价于书附录 C.1 的
    circ(sqrt(X^2+Y^2)/R)(Goodman 半径约定;表 1.2 印刷定义差 2 倍,
    见 ch1 勘误与第 6 章坑与陷阱)。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
from base import coords, rect, circle                     # noqa: E402  公共积木


# ---------- 频率坐标(书 coh_image 第 24 行) ----------
def freq_coords(L, M):
    """书第 24 行 fu=-1/(2*du):1/L:1/(2*du)-(1/L) 的等价实现(居中排列)。

    与 base.coords 的频率轴完全相同:间隔 1/L,范围 [-1/(2du), 1/(2du))。
    """
    return coords(L, M)[1]


# ---------- 相干传递函数(式 7.7:H = P(-zXP*lam*fU, -zXP*lam*fV)) ----------
def H_circ(fu, f0):
    """圆瞳相干传递函数(式 7.13):H = circ(sqrt(fU^2+fV^2)/f0)。

    f0 = wXP/(lam*zXP) 为相干截止频率(式 7.10);圆对称,负号可忽略。
    """
    FU, FV = np.meshgrid(fu, fu)
    return circle(FU, FV, f0)          # = circ(ρ/f0),Goodman 半径约定


def H_square(fu, f0):
    """方瞳相干传递函数(式 7.9):H = rect(fU/2f0)rect(fV/2f0)。

    式 7.9 的宗量 -lam*zXP*f/(2wXP) = -f/(2f0);rect 偶对称,负号可忽略。
    """
    FU, FV = np.meshgrid(fu, fu)
    return rect(FU / (2 * f0)) * rect(FV / (2 * f0))


def H_rect(fu, fv, f0x, f0y):
    """矩形瞳(习题 7.3):半宽 wx、wy 分别给出截止 f0x、f0y(式 7.10 推广)。"""
    FU, FV = np.meshgrid(fu, fv)
    return rect(FU / (2 * f0x)) * rect(FV / (2 * f0y))


def H_annulus(fu, f0o, f0i):
    """环瞳(习题 7.4):外半径 wo、内半径 wi 的环形区域,circ(r/wo)-circ(r/wi)。"""
    FU, FV = np.meshgrid(fu, fu)
    rho = np.hypot(FU, FV)
    return np.where(rho < f0o, 1.0, 0.0) - np.where(rho < f0i, 1.0, 0.0)


# ---------- 成像管线 ----------
def coh_image(ug, H):
    """相干成像(式 7.16,书 coh_image 第 35-39 行):
    ui = F^-1{ H * F{ug} }。返回 ui 与谱 Gi(居中排列)。"""
    Hf = np.fft.fftshift(H)                        # 书 35 行:移到 FFT 排列
    Gg = np.fft.fft2(np.fft.fftshift(ug))          # 书 36 行
    Gi = Gg * Hf                                   # 书 37 行
    ui = np.fft.ifftshift(np.fft.ifft2(Gi))        # 书 38 行
    return ui, np.fft.fftshift(Gi)                 # Gi 转回居中排列方便标注


def otf_from_H(H):
    """OTF 的 FFT 自相关算法(书 incoh_image 第 22-23 行,式 7.25/7.26)。

    OTF = ifft2(|fft2(fftshift(H))|^2),再除以零频值归一。
    注意:输出故意留在 FFT 排列(书第 22 行不 fftshift)——正好与
    fft2(fftshift(Ig)) 的谱排列一致,直接逐点相乘。
    """
    OTF = np.fft.ifft2(np.abs(np.fft.fft2(np.fft.fftshift(H))) ** 2)  # 书 22 行
    OTF = OTF / OTF[0, 0]                          # 书 23 行:除 OTF(1,1) 归一
    return OTF


def incoh_image(Ig, OTF):
    """非相干成像(式 7.33,书 incoh_image 第 32-36 行):
    Ii = F^-1{ OTF * F{Ig} }。
    数值上可能出小的虚部/负值(书第 35-36 行),取实部并把负值清零。"""
    Gg = np.fft.fft2(np.fft.fftshift(Ig))          # 书 32 行
    Gi = Gg * OTF                                  # 书 33 行(OTF 为 FFT 排列)
    Ii = np.fft.ifftshift(np.fft.ifft2(Gi))        # 书 34 行
    Ii = np.real(Ii)                               # 书 36 行:去残余虚部
    Ii = np.where(Ii >= 0, Ii, 0.0)                # 书 36 行:负值清零
    return Ii


# ---------- OTF 解析式 ----------
def otf_circ_analytic(fu, f0):
    """圆瞳衍射限 OTF(式 7.31,Goodman 图解自相关结果):
    H = (2/pi)[arccos(r) - r*sqrt(1-r^2)], r = rho/(2f0), rho <= 2f0;否则 0。"""
    FU, FV = np.meshgrid(fu, fu)
    rho = np.hypot(FU, FV)
    r = np.clip(rho / (2 * f0), 0.0, 1.0)
    H = (2 / np.pi) * (np.arccos(r) - r * np.sqrt(1.0 - r ** 2))
    return np.where(rho <= 2 * f0, H, 0.0)


def otf_square_analytic(fu, f0):
    """方瞳衍射限 OTF(式 7.29):H = Λ(fU/2f0)Λ(fV/2f0)。

    式 7.29 的宗量 lam*zXP*f/(2wXP) = f/(2f0);非相干截止 2f0(式 7.30)。"""
    FU, FV = np.meshgrid(fu, fu)
    tri_u = np.where(np.abs(FU) <= 2 * f0, 1 - np.abs(FU) / (2 * f0), 0.0)
    tri_v = np.where(np.abs(FV) <= 2 * f0, 1 - np.abs(FV) / (2 * f0), 0.0)
    return tri_u * tri_v


# ---------- 附录 C.6:udelta ----------
def udelta(x):
    """附录 C.6 单位样本 delta:x 取整数值(六位舍入后)处置 1。

    与 ucomb 同一家族:先把输入舍入到 10^-6 截断舍入误差,再判 ==0。
    习题 7.10 里点源 = udelta((X-x0)/dx) * udelta((Y-y0)/dy)。
    """
    x = np.round(np.asarray(x, dtype=float) * 1e6) / 1e6   # 书:round to 10^6ths
    return np.where(x == 0, 1.0, 0.0)                      # 书:out = x==0
