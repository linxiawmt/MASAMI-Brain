# -*- coding: utf-8 -*-
"""fourier 公共基础模块（Voelz《Computational Fourier Optics》笔记配套）。

这里集中了全书反复使用的基本函数与网格/变换约定：
  - 基本函数（书表 1.2）：rect / tri / sinc / gauss / circle
  - 圆盘频谱 jinc（书 1.5 节，注意半径约定，见 ch1 勘误 1）
  - 居中 DFT 定标包装：ft1_centered / ft2_centered（书 2.3 节）
  - 采样网格坐标 coords（书 2.4 节式 2.21 / 2.22）

来源说明：rect/tri/sinc/gauss/circle/jinc 与 ft*_centered 最初以内嵌函数
形式写在 python/fourier/ch1/ 的验证脚本里，ch2 起上移到本模块供后续章节
复用（定义与 ch1 页面内代码一致）。ch1 的脚本保持自包含、不改写。

import 路径处理（ch2/ch3... 各章脚本统一写法）：
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "common"))
    from base import rect, coords, ft2_centered

MATLAB/Python 差异（详见笔记第 2 章坑与陷阱）：
  - numpy 索引 0-based（MATLAB 1-based）：居中排列的"零坐标"样本，
    MATLAB 位置是 (N/2+1, M/2+1)，numpy 是 (N//2, M//2)。
  - np.meshgrid 默认 indexing='xy'，与 MATLAB meshgrid 行为一致；
    若在 Python 里改用 'ij'，两个维度会互换（圆对称函数无影响）。
"""

import numpy as np

__all__ = ["rect", "tri", "sinc", "gauss", "circle", "jinc",
           "ft1_centered", "ft2_centered", "coords"]


# ---------- 一维基本函数（书表 1.2） ----------

def rect(x):
    """rect(x)：|x| < 1/2 取 1，|x| = 1/2 取 1/2，其余取 0（表 1.2）。"""
    x = np.asarray(x, dtype=float)
    r = np.where(np.abs(x) < 0.5, 1.0, 0.0)
    return np.where(np.abs(x) == 0.5, 0.5, r)   # 边界半值；均匀网格很少命中


def tri(x):
    """Λ(x)：|x| <= 1 取 1-|x|，其余取 0（表 1.2）。"""
    x = np.asarray(x, dtype=float)
    return np.where(np.abs(x) <= 1.0, 1.0 - np.abs(x), 0.0)


def sinc(x):
    """归一化 sinc：sin(pi*x)/(pi*x)，x=0 处极限 1。
    与 np.sinc、MATLAB Signal Processing Toolbox 的 sinc 定义一致。"""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    nz = x != 0
    out[nz] = np.sin(np.pi * x[nz]) / (np.pi * x[nz])
    return out


def gauss(x):
    """一维高斯 exp(-pi*x^2)：傅里叶变换恰好是自己（表 1.3）。"""
    return np.exp(-np.pi * np.asarray(x, dtype=float) ** 2)


def circle(X, Y, R=0.5):
    """circ(sqrt(X^2+Y^2)/2R)：半径 R 的单位高度圆盘。
    表 1.2 的 circ(sqrt(x^2+y^2)) 对应 R = 1/2。
    注意表 1.3 的 circ 变换对按"半径 = 分母"计算（Goodman 约定），
    两种约定差 2 倍半径——数值验证见 ch1_verify_circ_jinc.py。"""
    return np.where(np.hypot(X, Y) < R, 1.0, 0.0)


def jinc(rho, R=1.0):
    """半径 R 圆盘的解析频谱：R*J1(2*pi*R*rho)/rho；rho=0 处极限 = 圆盘面积 pi*R^2。
    （jinc 的振幅版即圆孔夫琅禾费衍射图样。）"""
    rho = np.asarray(rho, dtype=float)
    from scipy.special import j1
    out = np.empty_like(rho)
    nz = rho != 0
    out[nz] = R * j1(2 * np.pi * R * rho[nz]) / rho[nz]
    out[~nz] = np.pi * R ** 2
    return out


# ---------- 网格坐标（书 2.4 节） ----------

def coords(L, N):
    """书式 (2.21)/(2.22) 的一维采样坐标。

    参数：L 采样空间边长 [m]，N 采样点数（书约定取偶数）。
    返回：
      x = n*dx,  n = -N/2 .. N/2-1,  dx = L/N      （式 2.21，居中升序）
      f = n/(N*dx) = n/L，即 -1/(2dx) .. 1/(2dx)-1/L（式 2.22，与 x 同序）

    返回的数组是"居中排列"（零坐标在中间）。喂给 fft 前，数据数组要先
    ifftshift 成"移位排列"（index 0 = 零坐标）；坐标轴本身只是标注用的
    升序序列，不参与移位——这是初学者最常见的搞混点，见笔记 2.4 节。
    """
    dx = L / N
    n = np.arange(N) - N // 2        # Python 0-based：n//2 即 MATLAB 的 N/2+1 号
    x = n * dx                       # x ∈ [-L/2, L/2 - dx]
    f = n / (N * dx)                 # f ∈ [-1/(2dx), 1/(2dx) - 1/L]
    return x, f


# ---------- 居中 DFT 定标包装（书 2.3 节） ----------

def ft1_centered(g, dx):
    """一维连续傅里叶变换的 DFT 近似（居中输出）：
        G(f) ≈ dx * Σ g(x) exp(-j2π f x)
    numpy.fft 的正变换不除 N（归一化在逆变换里），乘 dx 回到连续量纲。
    输入 g 按居中排列；内部 ifftshift -> fft -> fftshift。"""
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(g))) * dx


def ft2_centered(g, dx):
    """二维连续傅里叶变换的 DFT 近似（居中输出），原理同 ft1_centered，
    定标因子是 dx^2。"""
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(g))) * dx ** 2
