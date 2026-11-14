# ===== 书 5.1 · 相干模式 vs 非相干模式记录（式 5.1-6b / 5.1-6c，例 5.1 式 5.1-9）=====
# 演示广义双瞳处理器的核心结论（Poon & Indebetouw 2003）：
#   - 非相干模式（开掩膜 M=1）：记录的是强度物体 |T|^2，纯相位结构完全不可见（式 5.1-6c）；
#   - 相干模式（针孔掩膜 M=delta 且 p1=delta）：记录的是振幅物体 T 本身（式 5.1-6b），
#     p2=1 时就是例 5.1 的"复菲涅耳波带片全息图"（式 5.1-9）。
# 物体设计：三缝（振幅物体）+ 纯相位圆盘（|T|=1，相位 1.5 rad）。
#
# 严格书约定：前向变换 +j 指数 -> numpy 用 ifft2*N^2 实现；逆变换 -j 指数 -> fft2/N^2。
# （MATLAB 同理：ifft2(...)*N^2 与 fft2(...)/N^2，见同目录 .m 文件。）
#
# MATLAB / Python 差异提示：
# 1) numpy 的 meshgrid 默认 'xy'；MATLAB 的 meshgrid 同为 [X,Y]=meshgrid(x,y)，
#    两者数组均是 (行=y, 列=x)，一致。
# 2) plt.imshow 加 origin='lower' 让 y 轴向上（MATLAB imagesc 默认 y 向上）。

import numpy as np
import matplotlib.pyplot as plt

# ----- 0. 网格与书协议的傅里叶封装 -----
N = 512
dx = 1.0                      # 任意单位（无量纲演示）
lam0 = 0.5
k0 = 2 * np.pi / lam0
z0 = 2000.0                   # 物体前表面到焦点 c 的距离（式 5.1-1b 的 z+z0 取 z=0）

x = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(x, x)                     # (行=y, 列=x)
KX = 2 * np.pi * np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(KX, KX)
K2 = KX**2 + KY**2

def fwd(u):                                   # 书前向变换（+j 指数）
    return np.fft.ifft2(u) * N**2

def inv(U):                                   # 书逆变换（-j 指数）
    return np.fft.fft2(U) / N**2

def propagate(u, z):                          # 菲涅耳传播（式 2.3-13 传递函数）
    return inv(fwd(u) * np.exp(-1j * k0 * z) * np.exp(1j * z * K2 / (2 * k0)))

# ----- 1. 物体 T(x,y;z=0)：振幅三缝 + 纯相位圆盘 -----
T = np.zeros((N, N), dtype=complex)
for xc in (-180.0, -140.0, -100.0):           # 三条振幅缝（宽 12）
    T[(np.abs(X - xc) <= 6) & (np.abs(Y) <= 120)] = 1.0
disc = (X - 70.0)**2 + Y**2 <= 60.0**2         # 纯相位圆盘：零吸收，相位延迟 1.5 rad
T[disc] = np.exp(1j * 1.5)

# ----- 2. 相干模式记录（式 5.1-6b + 例 5.1 的 p2=1 -> 式 5.1-9）-----
# i_Omegap(x,y) = exp[-j k0 (z+z0)] * j k0 / (2 pi (z+z0)) * [ chirp(z+z0) * T ]
# 其中 chirp(Z) = exp[-j k0 (x^2+y^2) / (2Z)]，其频谱是 exp[+j Z k^2 / (2 k0)]。
A_coh = 1j * k0 / (2 * np.pi * z0) * np.exp(-1j * k0 * z0)          # 式 (5.1-9) 前因子
Coh = A_coh * inv(fwd(T) * np.exp(1j * z0 * K2 / (2 * k0)))          # 复 FZP 全息图

# 相干重建：把全息图当透射片，反向共轭传播即得 T（乘回前因子）
Rc = inv(fwd(Coh) * np.exp(-1j * z0 * K2 / (2 * k0))) / A_coh

# ----- 3. 非相干模式记录（式 5.1-6c，M=1；p1=1 球面波 + p2=delta 平面波，同 ch3 标准 OSH）-----
# i_Omegap = P1* P2 与 |T|^2 的相关 = |T|^2 与复共轭核的卷积；取 ch3 (3.5-5) 的 +chirp 分支：
# H = |T|^2 * exp[+j k0 r^2 / (2 z0)] * k0/(2 pi z0)，其频谱是 exp[-j z0 k^2 / (2 k0)]。
I2 = np.abs(T)**2
A_inc = k0 / (2 * np.pi * z0)
Inc = A_inc * inv(fwd(I2) * np.exp(-1j * z0 * K2 / (2 * k0)))        # 复 FZP 全息图（强度物体）

Ri = inv(fwd(Inc) * np.exp(1j * z0 * K2 / (2 * k0))) / A_inc         # 重建回 |T|^2

# ----- 4. 数值验证 -----
def corr(a, b):
    a = a.ravel() - a.mean(); b = b.ravel() - b.mean()
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))

ph_coh = np.angle(np.mean(Rc[disc]))            # 相干重建在圆盘上读出的相位
ph_inc = np.angle(np.mean(Ri[disc]))            # 非相干重建在圆盘上读出的相位
print('--- 式 5.1-6b（相干） vs 5.1-6c（非相干）记录对比 ---')
print(f'相干重建：corr(|Rc|, |T|)   = {corr(np.abs(Rc), np.abs(T)):+.4f}')
print(f'相干重建：圆盘相位读出     = {ph_coh:+.4f} rad（真值 1.5，纯相位结构被完整记录）')
print(f'非相干重建：corr(|Ri|,|T|^2) = {corr(np.abs(Ri), I2):+.4f}')
print(f'非相干重建：圆盘相位读出   = {ph_inc:+.4f} rad（只记录了 |T|^2，相位信息丢失）')
print('结论：掩膜从"开"换到"针孔"，记录对象从 |T|^2 变成 T —— 这就是 3-D 复相干处理的入口')

# ----- 5. 可视化（振幅为零处相位无定义，遮成黑色避免"angle(0±eps) 噪声"误导）-----
def phase_masked(U, floor=1e-9):
    ph = np.angle(U)
    ph[np.abs(U) < floor] = np.nan
    return ph

def arg_object(T):
    ph = np.angle(T)
    ph[np.abs(T) < 1e-9] = np.nan
    return ph

plt.rcParams['figure.facecolor'] = 'white'
fig, ax = plt.subplots(2, 3, figsize=(15, 9.2))
ext = [x[0], x[-1], x[0], x[-1]]
cmap = plt.get_cmap('jet').with_extremes(bad='k')
im = ax[0, 0].imshow(arg_object(T), extent=ext, cmap=cmap, vmin=-np.pi, vmax=np.pi)
ax[0, 0].set_title('1. Object arg(T): slits (amp) + disc (1.5 rad, |T|=1)')
plt.colorbar(im, ax=ax[0, 0], fraction=0.046)
ax[0, 1].imshow(np.abs(Coh), extent=ext, cmap='gray')
ax[0, 1].set_title('2. Coherent hologram |i_Omegap| (Eq.5.1-9, complex FZP)')
im = ax[0, 2].imshow(phase_masked(Rc), extent=ext, cmap=cmap, vmin=-np.pi, vmax=np.pi)
ax[0, 2].set_title('3. Coherent recon: arg(Rc) reads 1.5 rad on disc')
plt.colorbar(im, ax=ax[0, 2], fraction=0.046)

ax[1, 0].imshow(I2, extent=ext, cmap='gray', vmin=0, vmax=1.2)
ax[1, 0].set_title('4. Intensity object |T|^2 (what M=1 mode sees)')
ax[1, 1].imshow(np.abs(Ri), extent=ext, cmap='gray')
ax[1, 1].set_title('5. Incoherent recon: |Ri| ~ |T|^2 (slits + flat disc)')
im = ax[1, 2].imshow(phase_masked(Ri), extent=ext, cmap=cmap, vmin=-np.pi, vmax=np.pi)
ax[1, 2].set_title('6. Incoherent recon: arg(Ri) ~ 0 (no phase info)')
plt.colorbar(im, ax=ax[1, 2], fraction=0.046)
for a in ax.flat:
    a.set_xlabel('x'); a.set_ylabel('y')
plt.tight_layout()
plt.savefig('sim-coherent-incoherent.png', dpi=150, bbox_inches='tight')
plt.show()
