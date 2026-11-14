# ===== 书 3.6 节 · Table 3.2 TDFZP.m 的 Python 版：时间依赖菲涅耳波带片的运行条纹 =====
# TDFZP(x,y;z0,t) = A + B*sin[(k0/2*z0)(x^2+y^2) + Omega*t]     (3.6-1)
# 平面波（w0+Omega）与球面波（w0）干涉形成的"会跑的 FZP"：t 增大时同相位圆环半径缩小，
# 条纹向中心运动（running fringes toward the center）——OSH 就是用它光栅扫描 3-D 物体。
# 换瞳（p1=delta, p2=1）得 TDFZP(x,y;z0,-t)：条纹远离中心；sin 全息图不变、cos 全息图反号 (3.6-4)。
#
# 原书 TDFZP.m 参数：B=10.01e6（时间频率，任意）、D=6（空间尺度，任意）、t=linspace(0,1,35)，
# x,y=linspace(-2.5,2.5,256)，FZP = 1+sin(D*(x^2+y^2)+B*t)，逐帧归一化后 movie 播放。
# 本版固定 3 帧（t 使相位差可见）输出静态对比；想看动画可仿照原书逐帧绘制。
#
# MATLAB / Python 差异提示：
# 1) 原书三重循环 (ii,jj,kk) 用矩阵广播一行完成：D*(X^2+Y^2)+B*t。
# 2) 原书 image(256*FZP) 行方向向下；plt.imshow 加 origin='lower' 对齐数学坐标。
# 3) 改变条纹运动方向：把 sin 里的 +B*t 改成 -B*t（等价于 3.6-4 的 TDFZP(x,y;z0,-t)）。

import numpy as np
import matplotlib.pyplot as plt

Dc = 6.0                       # 空间尺度（原书 D=6，任意单位）
Bf = 10.01e6*2*np.pi*0.05      # 时间相位速度（原书 B=10.01e6，此处缩放使静态帧相位差可见）
xg = np.linspace(-2.5, 2.5, 256)
Xg, Yg = np.meshgrid(xg, xg)

fig, ax = plt.subplots(1, 3, figsize=(12, 4))
for i, tt in enumerate([0.0, 0.35, 0.7]):
    FZP = 1 + np.sin(Dc*(Xg**2 + Yg**2) + Bf*tt)     # TDFZP (3.6-1)
    ax[i].imshow(FZP, cmap='gray', extent=[-2.5, 2.5, -2.5, 2.5], origin='lower')
    ax[i].set_title(f'TDFZP: sin[6r$^2$ + $\\Omega$t],  $\\Omega$t = {Bf*tt:.2f}')
    ax[i].axis('off')
plt.tight_layout(); plt.show()

# 条纹向中心运动的验证：跟踪相位 = D*r^2 + B*t = const 的半径 r(t) = sqrt((c - B*t)/D)，
# t 增大 -> r 减小 -> 圆环向中心收缩（书 3.6 节 "running fringes moving toward the center"）。
