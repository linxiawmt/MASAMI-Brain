# -*- coding: utf-8 -*-
"""第 5 章三个传播器:propTF / propIR / propFF(书 5.1/5.2/5.5 节)。

与 matlab/fourier/ch5/propTF.m / propIR.m / propFF.m 逐行对应,供本章
各脚本(ch5_sqr_beam / ch5_sampling / ch5_fraunhofer / ch5_exercises)导入。
固定舞步:"居中排列 <-> FFT 排列"之间靠 fftshift 进、ifftshift 回
(偶数 M 时 numpy 的 fftshift 与 ifftshift 相同,书全用偶数 M)。

MATLAB/Python 差异(详见笔记第 5 章坑与陷阱):
  - numpy 索引 0-based(MATLAB 1-based):居中排列的"零坐标"样本,
    MATLAB 位置是 (M/2+1, M/2+1),numpy 是 (M//2, M//2)。
  - numpy.fft 与 MATLAB fft 定标一致(正变换不除 N、逆变换除 N),
    propIR 里 U1 的 FFT 与 U2 的 IFFT 定标互相抵消,唯独 FFT{h} 需要
    自己的 dx^2 才能回到连续量纲(书 5.2 节注 b)。
  - propTF/propIR 的源面与观察面共边长(L1 = L2 = L);propFF 的观察面
    边长由源面参数决定(式 5.18),Python 版返回 (u2, L2) 元组。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
from base import coords                              # noqa: E402  公共积木(ch1 起)


# ===== 书 5.1 节:propTF(与 matlab/fourier/ch5/propTF.m 逐行对应) =====
def propTF(u1, L, lam, z):
    """propagation - transfer function approach(书 5.1 节,式 5.1/5.2)。

    假设源面与观察面边长相同(L1 = L2 = L)、均匀采样。
    H = exp(jkz)*exp[-j pi lam z (fX^2+fY^2)] 中的整体相位 exp(jkz)
    不影响横向结构,按书第 18 行的做法直接丢掉。
    """
    M, N = u1.shape                                   # 书 11 行:取阵列尺寸(只用 M)
    dx = L / M                                        # 书 12 行:采样间隔
    k = 2 * np.pi / lam                               # 书 13 行:波数
    _, fx = coords(L, M)                              # 书 15 行:频率坐标(居中)
    FX, FY = np.meshgrid(fx, fx)                      # 书 16 行(fy 与 fx 相同,省一行)
    H = np.exp(-1j * np.pi * lam * z * (FX ** 2 + FY ** 2))   # 书 18 行:传递函数(丢 e^{jkz})
    H = np.fft.fftshift(H)                            # 书 19 行:移位到 FFT 排列
    U1 = np.fft.fft2(np.fft.fftshift(u1))             # 书 20 行:源场移位 + FFT
    U2 = H * U1                                       # 书 21 行:逐点相乘
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))           # 书 22 行:逆 FFT + 移回居中
    return u2


# ===== 书 5.2 节:propIR(与 matlab/fourier/ch5/propIR.m 逐行对应) =====
def propIR(u1, L, lam, z):
    """propagation - impulse response approach(书 5.2 节,式 5.3/5.4)。

    与 propTF 解析上等价,但离散采样后行为不同(书 5.4 节的主题):
      - h 的 chirp 相位欠采样(Dx > lam z/L)时,FFT 的周期延拓会在
        观察面产生周期副本;
      - h 过采样(Dx < lam z/L)时,FFT 本身起到"窗口滤波"作用,
        压掉 L/(2 lam z) 以外的源频谱——长距离下反而更平滑。
    """
    M, N = u1.shape                                   # 书 11 行
    dx = L / M                                        # 书 12 行
    k = 2 * np.pi / lam                               # 书 13 行
    x, _ = coords(L, M)                               # 书 15 行:空间坐标(居中)
    X, Y = np.meshgrid(x, x)                          # 书 16 行
    h = 1 / (1j * lam * z) * np.exp(1j * k / (2 * z) * (X ** 2 + Y ** 2))  # 书 18 行:脉冲响应
    H = np.fft.fft2(np.fft.fftshift(h)) * dx ** 2     # 书 19 行:FFT 得传递函数,
    #                                                #   乘 dx^2 回到连续量纲——
    #                                                #   U1 与 U2 的 FFT/IFFT 定标互相抵消,
    #                                                #   只有 FFT{h} 这一步需要自己的定标(书 5.2 节注 b)
    U1 = np.fft.fft2(np.fft.fftshift(u1))             # 书 20 行
    U2 = H * U1                                       # 书 21 行
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))           # 书 22 行
    return u2


# ===== 书 5.5 节:propFF(与 matlab/fourier/ch5/propFF.m 逐行对应) =====
def propFF(u1, L1, lam, z):
    """propagation - Fraunhofer pattern(书 5.5 节,式 5.16-5.19)。

    观察面边长与采样间隔由源面参数决定(式 5.18):
        L2  = lam*z/dx1,   dx2 = lam*z/L1
    所以函数额外返回 L2,省得主脚本再算一遍(书原话)。
    """
    M, N = u1.shape                                   # 书 11 行
    dx1 = L1 / M                                      # 书 12 行:源面采样间隔
    k = 2 * np.pi / lam                               # 书 13 行
    L2 = lam * z / dx1                                # 书 15 行:观察面边长
    dx2 = lam * z / L1                                # 书 16 行:观察面采样间隔
    x2, _ = coords(L2, M)                             # 书 17 行:观察面坐标(居中,与书式 5.19 一致)
    X2, Y2 = np.meshgrid(x2, x2)                      # 书 18 行
    c = 1 / (1j * lam * z) * np.exp(1j * k / (2 * z) * (X2 ** 2 + Y2 ** 2))  # 书 20 行:前导因子
    u2 = c * np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u1))) * dx1 ** 2   # 书 21 行
    return u2, L2
