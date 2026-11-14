# ===== 书 2.3.2 · Table 2.1 Fresnel_diffraction.m 的 Python 版 =====
# 正方形孔径的菲涅耳衍射仿真（卷积路线，式 (2.3-12)）：
#     psi_p(x,y;z) = psi_p0(x,y) * h(x,y;z)
# 输入：M x M 正方形孔 rect(x/a, y/a)，a = M/256 = 0.4336 cm
# 光源：He-Ne 红光 632.8 nm，观察屏 z = 5 cm（复现原书图 2.11）
#
# 与第 1 章拓展 1.2.3 的 fresnel_lsi.py 是同一物理的"书中正式版"：
# 那里把它当"自由空间 LSI 系统"的练习，这里按书 Table 2.1 原流程实现。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 矩阵 1-based、Python 0-based：MATLAB 第 127/128 列 = Python 下标 126/127。
# 2) 本脚本用"显式构造 h 再做频域乘法"的卷积路线——卷积定理只依赖同一对
#    正/逆变换，与 +j/-j 指数约定无关，所以照抄书中 fft2 写法结果就是书的卷积
#    （h 本身已按书的 (2.3-11) 用 -j 二次相位构造）。"约定坑"只在用解析
#    传递函数 H(kx,ky;z)（式 2.3-13）时才咬人，见 asm_propagation.py。
# 3) exp(-4*200*r^2) 是高斯窗：把脉冲响应截断在显示窗口内，抑制离散卷积的
#    卷绕误差（原书写法，4*200=800）。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 参数（与原书 Table 2.1 完全一致） -----
L = 1.0                    # L : length of display area（显示区域 L x L，单位 cm）
N = 256                    # N : number of sampling points（采样点数）
dx = L / (N - 1)           # dx : step size（采样步长 1/255 cm）
M = 111                    # 正方形孔的像素宽度（a = M/256 = 0.4336 cm）
lam = 0.6328e-4            # He-Ne 红光波长（cm）
z = 5.0                    # 观察距离（cm）

# ----- 2. 输入：中心 M x M 正方形孔 rect(x/a, y/a) -----
R = np.zeros((N, N))       # assign a matrix (256x256) of zeros
n = (M - 1) // 2
c = 128 - 1                # MATLAB 的第 128 个元素（1-based）= Python 下标 127
R[c - n:c + n + 1, c - n:c + n + 1] = 1.0   # 背景置 0、中心置 1 的方孔

# ----- 3. 坐标轴标定（与书中 fft2Drect.m 相同的方式） -----
k = np.arange(N)
X = k / (N - 1) - L / 2    # 空间坐标 [-0.5, 0.5] cm
Y = X.copy()

# ----- 4. 菲涅耳脉冲响应 h(x,y;z)（式 2.3-11 + 高斯窗） -----
sigma = np.pi / (lam * z)                  # sigma = pi/(lambda*z)，二次相位系数
r2 = X[:, None]**2 + Y[None, :]**2         # x^2+y^2 网格（广播代替双层 for）
h = 1j * (sigma / np.pi) * np.exp(-1j * sigma * r2) * np.exp(-800.0 * r2)
#   j*(sigma/pi) = j/(lambda*z) 振幅因子（全局相位 exp(-j k0 z) 未算，只影响整体相位）

# ----- 5. 频域实现卷积（照抄原书流程：FR.*H 后 ifft2 再 fftshift） -----
FR = np.fft.fft2(R) / N**2                 # 原书：FR=(1/256)^2*fft2(R)
H = np.fft.fft2(h) / N**2                  # 原书：H=(1/256)^2*fft2(h)
hr = np.fft.ifft2(FR * H) * N**2           # 原书：hr=ifft2(HR); hr=(256^2)*hr
hr = np.fft.fftshift(hr)                   # 把卷积中心搬回图像中央
I = np.abs(hr)**2                          # 光强

# ----- 6. 可视化（对应原书 figure(1)(2)(3) 与 CCD 视图） -----
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))

ax[0, 0].imshow(R, extent=[X[0], X[-1], Y[0], Y[-1]], cmap='gray', origin='lower')
ax[0, 0].set_title('Input square aperture rect(x/a, y/a)')
ax[0, 0].set_xlabel('x (cm)'); ax[0, 0].set_ylabel('y (cm)')

ax[0, 1].plot(X + dx / 2, R[:, 126], '-o', ms=2)      # MATLAB R(:,127) -> Python 第 126 列
ax[0, 1].set_title('Square aperture: cross-section')
ax[0, 1].set_xlabel('x (cm)'); ax[0, 1].axis([-0.5, 0.5, -0.1, 1.2]); ax[0, 1].grid(True)

ax[1, 0].plot(X + dx / 2, np.abs(hr[:, 126]))         # 衍射幅度中心横截面（原图 2.11c）
ax[1, 0].set_title('|psi_p(x,y;z=5cm)| cross-section')
ax[1, 0].set_xlabel('x (cm)'); ax[1, 0].grid(True)

ax[1, 1].imshow(I, extent=[X[0], X[-1], Y[0], Y[-1]], cmap='gray', origin='lower')
ax[1, 1].set_title('2-D diffraction pattern (CCD view) at z = 5 cm')
ax[1, 1].set_xlabel('x (cm)'); ax[1, 1].set_ylabel('y (cm)')

plt.tight_layout()
plt.show()
