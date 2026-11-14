# ===== 拓展 1.2.3 · 自由空间当作 LSI 系统：菲涅耳衍射的频域实现 =====
# 改编自你已有的 Fresnel Diffraction.m（规范化重写，物理参数与结果一致）：
# 矩形孔 f(x,y) 放在 z=0，求 z=5 cm 处的衍射光强。
# 本章视角（1.2 节）：自由空间是一个 LSI 系统，
#     路线 A（空间域卷积）：g = f * h，h 为菲涅耳脉冲响应；
#     路线 B（频域乘传递函数）：G = F*H —— 本脚本用的就是卷积定理的频域实现。
#
# MATLAB / Python 差异提示：
# 1) 离散卷积定理实现里，FFT 归一化因子 (1/N^2) 配 (N^2) 成对出现，互相抵消；
#    MATLAB 与 numpy 的 fft2/ifft2 约定一致，可逐行照抄。
# 2) meshgrid 默认 indexing='xy'：h(r,c) 中 MATLAB 的 X(r),Y(c) 对应
#    Python 里先对 X、Y 各取一列/一行组合，这里直接用广播 X[:,None]**2+Y[None,:]**2。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 参数（与 Fresnel Diffraction.m 完全一致） -----
L = 1.0                    # 显示区域 L x L（单位 cm）
N = 256                    # 采样点数
dx = L / (N - 1)           # 采样步长 1/255 cm
M = 111                    # 正方形孔的像素宽度（111*dx ≈ 0.43 cm）
lam = 6.328e-5             # He-Ne 红光波长（cm），即 632.8 nm
z = 5.0                    # 传播距离（cm）

# ----- 2. 输入：中心 M x M 正方形孔 -----
R = np.zeros((N, N))
n = (M - 1) // 2
c = 128 - 1                # MATLAB 第 128 个元素 = Python 下标 127
R[c - n:c + n + 1, c - n:c + n + 1] = 1.0

# ----- 3. 坐标轴（与书中 fft2Drect.m 相同的标定方式） -----
k = np.arange(N)
X = k / (N - 1) - L / 2            # 单位 cm
Y = X.copy()

# ----- 4. 菲涅耳脉冲响应 h(x,y;z)（二次相位 + 高斯窗防混叠） -----
sigma = np.pi / (lam * z)                  # σ = pi/(lambda*z)，二次相位系数
r2 = X[:, None]**2 + Y[None, :]**2         # x^2+y^2 网格
h = 1j * (sigma / np.pi) * np.exp(-1j * sigma * r2) * np.exp(-800.0 * r2)
#   ↑ j*(σ/π)=j/(λz) 振幅因子；exp(-800*r2) 是高斯低通窗（原书 4*200=800），
#     用来把脉冲响应截断在显示窗口内，防止离散卷积的卷绕误差。

# ----- 5. 频域实现卷积（卷积定理：F{f*h} = F*H，照抄原书 MATLAB 流程） -----
FR = np.fft.fft2(R) / N**2          # 原书：FR=(1/256)^2*fft2(R)
H = np.fft.fft2(h) / N**2           # 原书：H=(1/256)^2*fft2(h)
g = np.fft.fftshift(np.fft.ifft2(FR * H)) * N**2
# 原书：hr=ifft2(HR); hr=(256^2)*hr; hr=fftshift(hr)
# 说明 1：本书前向变换为"正指数"，严格写法是 ifft2*N^2；卷积定理下两条路线
#         只差共轭/翻转，对本例的中心对称输入结果相同（数值验证见 verify_ch1.py）。
# 说明 2：这是"自然顺序"的循环卷积，奇数宽孔径中心在 MATLAB 第 128 个元素，
#         卷积峰会绕回后落在中心偏左约 2 像素处——原书代码同样如此，显示无碍。

I = np.abs(g)**2                            # 光强

# ----- 6. 可视化（对应原书 figure(1)(2)(3)(5)） -----
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
ax[0, 0].imshow(R, extent=[X[0], X[-1], Y[0], Y[-1]], cmap='gray', origin='lower')
ax[0, 0].set_title('Input square aperture rect(x/a, y/a)')
ax[0, 0].set_xlabel('x (cm)'); ax[0, 0].set_ylabel('y (cm)')

ax[0, 1].plot(X + dx / 2, R[:, 126], '-o', ms=2)
ax[0, 1].set_title('rect function: cross-section plot')
ax[0, 1].set_xlabel('x (cm)'); ax[0, 1].axis([-0.5, 0.5, -0.1, 1.2]); ax[0, 1].grid(True)

ax[1, 0].plot(X + dx / 2, np.abs(g[:, 126]))
ax[1, 0].set_title('|g(x,y)| cross-section at z = 5 cm')
ax[1, 0].set_xlabel('x (cm)'); ax[1, 0].grid(True)

ax[1, 1].imshow(I, extent=[X[0], X[-1], Y[0], Y[-1]], cmap='gray', origin='lower')
ax[1, 1].set_title('2-D diffraction pattern (CCD view) at z = 5 cm')
ax[1, 1].set_xlabel('x (cm)'); ax[1, 1].set_ylabel('y (cm)')

plt.tight_layout()
plt.show()
