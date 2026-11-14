# ===== 书 4.1 节 · 分辨率-景深度量衡（式 4.1-1 ~ 4.1-5）与 Schilling 荧光实验的数值账 =====
# 量子力学不确定关系（书用 h 而非 hbar）：
#   横向：Δr·Δp_r >= h,  Δp_r = p0·sin(θ/2),  p0 = h/λ0
#       -> Δr >= λ0/sin(θ/2) = λ/(n0 sin(θ/2)) = λ/NA          (4.1-2)
#   纵向：Δz·Δp_z >= h,  Δp_z = p0·[1-cos(θ/2)]
#       -> Δz >= λ0/[1-cos(θ/2)] ≈ 2λ0/sin²(θ/2) = 2n0λ/NA²    (4.1-4)
#   两式结合：(Δr)²/Δz >= λ0/2                                  (4.1-5，标度律)
# 分辨率提高 2 倍 -> 景深必须缩短 4 倍。
#
# 实验参数（书 4.1 节，Schilling et al. Optics Lett. 22, 1506 (1997)）：
#   f = 150 mm（聚焦透镜 L1 焦距），D = 10 mm（被聚焦的平面波直径），λ0 = 514.5 nm（氩离子）
# 书印数值：NA ≈ 0.033，Δr ≈ 18.5 µm，Δz ≈ 1028.4 µm。
#
# MATLAB / Python 差异提示：纯数值脚本，无傅里叶约定问题。

import sys
import numpy as np
import matplotlib.pyplot as plt

if hasattr(sys.stdout, 'reconfigure'):      # Windows 控制台中文/µm 输出
    sys.stdout.reconfigure(encoding='utf-8')

# ----- 1. 实验的 NA 与衍射极限（式 4.1-2 / 4.1-4） -----
f1 = 150e-3          # L1 焦距 (m)
D = 10e-3            # 平面波直径 (m)
lam = 514.5e-9       # 氩离子 514.5 nm 线（真空中）
n0 = 1.0             # 样本空间近似空气/水，取 1

half_angle = np.arctan((D / 2) / f1)     # θ/2
NA = n0 * np.sin(half_angle)
d_r = lam / NA                           # (4.1-2) 横向分辨率
d_z = 2 * n0 * lam / NA**2               # (4.1-4) 景深

print(f'θ/2 = {np.degrees(half_angle):.4f} deg,  NA = {NA:.4f}')
print(f'式(4.1-2): Δr = λ0/NA = {d_r*1e6:.1f} µm   （书印 18.5 µm）')
print(f'式(4.1-4): Δz = 2λ0/NA² = {d_z*1e6:.1f} µm = {d_z*1e3:.2f} mm（书印 1028.4 µm）')
print(f'15 µm 荧光小球 vs Δr = {d_r*1e6:.1f} µm：接近横向极限，可分辨')
print(f'式(4.1-5) 标度律：(Δr)²/Δz = {(d_r**2/d_z)*1e6:.4f} µm ≈ λ0/2 = {lam/2*1e6:.4f} µm')

# ----- 2. (4.1-5) 标度律曲线：Δz(Δr) 双曲线 -----
r = np.linspace(5e-6, 60e-6, 400)
z_hyp = 2 * r**2 / lam                   # 由 (Δr)²/Δz = λ0/2 反解

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(r * 1e6, z_hyp * 1e3, 'b-', lw=1.5, label=r'$(\Delta r)^2/\Delta z=\lambda_0/2$  (4.1-5)')
ax[0].plot(d_r * 1e6, d_z * 1e3, 'ro', ms=8, label=f'Schilling setup ({d_r*1e6:.1f} µm, {d_z*1e3:.2f} mm)')
ax[0].annotate('resolution x2\n=> DOF /4', xy=(d_r * 1e6 / 2, d_z * 1e3 / 4),
               xytext=(15, 3), fontsize=9, color='r',
               arrowprops=dict(arrowstyle='->', color='r'))
ax[0].set_xlabel(r'lateral resolution $\Delta r$ (µm)')
ax[0].set_ylabel(r'depth of focus $\Delta z$ (mm)')
ax[0].set_title('Uncertainty trade-off in microscopy (Eq. 4.1-5)')
ax[0].grid(alpha=0.3); ax[0].legend(fontsize=9)

# ----- 3. FZP 条纹间距账：TDFZP 扫描场最细条纹（决定采样要求） -----
x = np.linspace(0, D / 2, 300)                     # 扫描场半径（受平面波直径限制）
fringe_period = 2 * np.pi / (2 * np.pi * x / (lam * f1))   # 局部周期 λ0·f1/x = λ0/(NA 的局域值)
ax[1].semilogy(x * 1e3, fringe_period * 1e6, 'g-', lw=1.5)
ax[1].axhline(2 * d_r * 1e6, color='r', ls='--', label=f'2·Δr = {2*d_r*1e6:.1f} µm')
ax[1].set_xlabel('radius in scan field (mm)')
ax[1].set_ylabel('local FZP fringe period (µm)')
ax[1].set_title('TDFZP fringe period across scan field')
ax[1].grid(alpha=0.3, which='both'); ax[1].legend(fontsize=9)
plt.tight_layout()
plt.savefig('sim-shm-resolution.png', dpi=150, bbox_inches='tight')
plt.show()

# 结论：NA=0.0333 时横向 Δr=15.5 µm、景深 Δz=0.93 mm——与书印 18.5 µm / 1028.4 µm
# 同数量级（书值对应 NA≈0.028 与 0.0316，属舍入不一致，见笔记拓展注）。
