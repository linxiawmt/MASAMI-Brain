# ===== 书 1.1 例 1.1 · Table 1.3 fft2Drect.m 的 Python 版 =====
# 二维 rect 函数（M x M 正方形）的傅里叶变换仿真，验证变换对
#     rect(x/a, y/b)  <=>  ab sinc(a*kx/2pi, b*ky/2pi)
# 书中运行时用 input 让你输 M（奇数 3~33），书中示例 M=11：
#     a = M/256 = 0.0429，第一零点 kx0 = 2*pi/a = 146.23 rad/(unit of length)
# 本版固定 M=11（想换值改这一个数即可）。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 矩阵 1-based、Python 0-based：MATLAB 的第 128 行 = Python 下标 127。
# 2) 书中前向变换按约定是"正指数"，严格实现应为 ifft2(R)*N^2；
#    但 rect 是实偶函数，fft2(R)/N^2 与之完全相等，所以书里直接用 fft2。
#    numpy.fft.fft2 与 MATLAB fft2 指数约定一致（负指数），这里照抄书里写法。
# 3) MATLAB image() 的行号方向与 plt.imshow 不同，用 origin='lower' 对齐数学坐标。

import numpy as np
import matplotlib.pyplot as plt

# ----- 显示区域与采样 -----
L = 1.0                    # display area is L by L, L has unit of length
N = 256                    # number of sampling points
dx = L / (N - 1)           # dx : step size

# ----- 生成 M x M 正方形 rect(x/a, y/a)，M=奇数 -----
M = 11                     # 书中示例值：input('M (size of rect(x/a,y/a), enter odd numbers from 3-33-)=')
a = M / 256
kx0 = 2 * np.pi / a
print(f'a= {a:0.5g} [unit of length]')                    # 对应书中 sprintf
print(f'kx0 (first zero)= {kx0:0.5g} [radian/unit of length]')

R = np.zeros((256, 256))   # assign a matrix (256x256) of zeros
n = (M - 1) // 2
c = 128 - 1                # MATLAB 的第 128 个元素（1-based）= Python 下标 127
R[c - n:c + n + 1, c - n:c + n + 1] = 1.0   # 背景置 0、中心置 1 的正方形孔

# ----- 坐标轴标定（书中的写法：零频移到中央） -----
k = np.arange(256)
X = (k / 255) * 1.0 - L / 2        # X(k)=1/255*(k-1)-L/2，即 (k-1)*dx-L/2
Y = X.copy()
# Kx(k)=2*pi*(k-1)/((N-1)*dx) - pi/dx；本例 (N-1)*dx=1
Kx = 2 * np.pi * k / ((N - 1) * dx) - (2 * np.pi * (256 - 1)) / ((N - 1) * dx) / 2
Ky = Kx.copy()

fig, ax = plt.subplots(2, 3, figsize=(13, 8))

# ----- 书 figure(1)：rect 函数灰度图 -----
im0 = ax[0, 0].imshow(R, extent=[X[0] + dx / 2, X[-1] + dx / 2, Y[0] + dx / 2, Y[-1] + dx / 2],
                      cmap='gray', origin='lower')
ax[0, 0].set_title('rect function: gray-scale plot')
ax[0, 0].set_xlabel('x'); ax[0, 0].set_ylabel('y')

# ----- 书 figure(2)：rect 的中心横截面 -----
ax[0, 1].plot(X + dx / 2, R[:, 126])       # MATLAB R(:,127) -> Python 第 126 列
ax[0, 1].set_title('rect function: cross-section plot')
ax[0, 1].set_xlabel('x'); ax[0, 1].set_ylabel('rect(x/a)')
ax[0, 1].grid(True)
ax[0, 1].axis([-0.5, 0.5, -0.1, 1.2])

# ----- 计算傅里叶变换（书：FR=(1/256)^2*fft2(R) 后 fftshift） -----
FR = np.fft.fft2(R) / 256**2
FR = np.fft.fftshift(FR)

# ----- 书 figure(3)：|F| 沿 kx 的横截面（幅度 x10 显示） -----
ax[0, 2].plot(Kx - np.pi / (dx * (N - 1)), 10 * np.abs(FR[:, 126]))
ax[0, 2].set_title('Square-absolute value of Fourier transform\nof rect function: cross-section plot')
ax[0, 2].set_xlabel('kx'); ax[0, 2].set_ylabel('|a*b*sinc(a*kx/2pi)|')
ax[0, 2].axis([-800, 800, 0, np.abs(FR).max() * 10.1])
ax[0, 2].grid(True)

# ----- 书 figure(4)：|F|^2 三维网格图（matplotlib 用等高线图代替 mesh） -----
ax[1, 0].contourf(Kx, Ky, np.abs(FR)**2, levels=50, cmap='viridis')
ax[1, 0].set_title('Square-absolute value of FT: 3-D plot (contour), scale arbitrary')
ax[1, 0].set_xlabel('kx'); ax[1, 0].set_ylabel('ky')
ax[1, 0].set_aspect('equal')

# ----- 书 figure(5)：|F|^2 灰度图（gain=10000 强压动态范围） -----
gain = 10000
im1 = ax[1, 1].imshow(np.abs(FR)**2 / np.abs(FR).max()**2,
                      extent=[Kx[0], Kx[-1], Ky[0], Ky[-1]],
                      cmap='gray', origin='lower', vmax=255.0 / gain)
ax[1, 1].set_title('Square-absolute value of FT of the rect\nfunction: gray-scale plot')
ax[1, 1].set_xlabel('kx'); ax[1, 1].set_ylabel('ky')

ax[1, 2].axis('off')       # 占位空面板
plt.tight_layout()
plt.show()
