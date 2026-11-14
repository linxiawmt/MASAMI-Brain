# ===== 书 2.5.1 · Table 2.2 Fresnel_zone_plate.m 的 Python 版 =====
# 画正弦菲涅耳波带片（FZP，点源全息图，式 (2.5-3)/(2.5-4)）：
#     FZP(x,y;z0) = A + B * sin( sigma * ((x-x0)^2 + (y-y0)^2) ),  sigma = 1/z0
# 复现原书图 2.19（z0 与 2z0 的深度依赖）与图 2.20（横向位置 (x0,y0) 依赖）。
# 所有尺度沿用原书"任意单位"（sigma=1/z0，x 步长 0.1，范围 ±12.8）；
# 显示 127*(1+imag(fFZP)) 即 127*(1+sin(sigma*r'^2))，灰度 0~254。
#
# MATLAB / Python 差异提示：
# 1) 原书用 input() 让你现场输入 z0（2~10）与 x0（-8~8），本版固定四组参数
#    一次性复现图 2.19 与 2.20；想改值改 z_list / x_off 即可。
# 2) 原书 image(R) 的行方向与 plt.imshow 默认相反，加 origin='lower' 对齐数学坐标。

import numpy as np
import matplotlib.pyplot as plt

ROWS = COLS = 256
s = np.arange(ROWS) * 0.1 - 12.8           # 原书 x=-12.8:0.1:12.7（步长 0.1，共 256 点）
X, Y = np.meshgrid(s, s)

def make_fzp(z0, x0, y0):
    """1+sin(sigma*((x-x0)^2+(y-y0)^2))，sigma=1/z0（原书任意单位）"""
    sigma = 1.0 / z0
    phase = sigma * ((X - x0)**2 + (Y - y0)**2)
    fFZP = np.exp(1j * phase)              # 原书：exp(j*sigma*(...)+j*sigma*(...))
    R = 127 * (1 + np.imag(fFZP))          # 归一化后取虚部 -> 1+sin(phase)
    return R

z1 = 4.0                                    # 基准深度（原书 input z0，2~10 内取值）
fig, ax = plt.subplots(2, 2, figsize=(10.5, 9.5))

# ----- 图 2.19 复现：同一点源（on-axis），深度 z0 与 2*z0 -----
im0 = ax[0, 0].imshow(make_fzp(z1, 0, 0), cmap='gray', origin='lower')
ax[0, 0].set_title(f'FZP, on-axis, z0 = {z1:g}  (lower fringe freq. far away)')

ax[0, 1].imshow(make_fzp(2 * z1, 0, 0), cmap='gray', origin='lower')
ax[0, 1].set_title(f'FZP, on-axis, z0 = {2*z1:g}  (deeper -> coarser fringes)')

# ----- 中心 1-D 剖面对比（原书图 2.19 下排 "1-D plots through center" -----
r1 = make_fzp(z1, 0, 0)
r2 = make_fzp(2 * z1, 0, 0)
ax[1, 0].plot(s, r1[:, COLS // 2] / 254.0, label=f'z0 = {z1:g}')
ax[1, 0].plot(s, r2[:, COLS // 2] / 254.0, label=f'z0 = {2*z1:g}')
ax[1, 0].set_title('1-D plots of FZP through center')
ax[1, 0].set_xlabel('x (a.u.)'); ax[1, 0].legend(loc='upper right')

# ----- 图 2.20 复现：同深度，横向位置从 (0,0) 挪到 (x0,y0) -----
x0 = 5.0
ax[1, 1].imshow(make_fzp(z1, x0, -x0), cmap='gray', origin='lower')   # 原书 y0=-x0
ax[1, 1].set_title(f'FZP center moves to (x0,y0)=({x0:g},{-x0:g}), same depth')

plt.tight_layout()
plt.show()
