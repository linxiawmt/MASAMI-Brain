# ===== 拓展实战（编者，原书无）· 同轴数字全息下 振幅物体 vs 相位物体 的成像对比 =====
# 改编自你已有的 Amplitude_and_Phase_Object.m / 相位物体和振幅物体.py（规范化重写，去 OpenCV 依赖）：
#   左：纯振幅粒子（强吸收 0.1，无相位）；右：纯相位粒子（无吸收，相位延迟 1.5 rad）
#   -> 角谱正向传播 5 cm 生成同轴全息图 -> 相机面取强度 -> sqrt 强度当复振幅估计 -> 反向传播重建。
#
# 与 4.1 节的关系（为什么扫描全息显微重要）：
# - 生物样品多为"相位物体"：几乎不吸收，只推迟相位——普通强度成像里近乎透明（本跑振幅通道
#   的相位粒子对比度只有 ~x%），重建振幅图里相位粒子只剩一圈"幽灵环"（孪生像干涉纹）。
# - 相位信息藏在同轴干涉条纹里：sqrt 强度 + 反向传播只能部分恢复（相位穹顶被背景干涉压低）。
# - SHM（4.1）记录的电子全息图是复值的（ic/is 双通道可合成复全息图，ch3 式 3.5-5），
#   相位通道一就读出定量相位——这正是 ch5"相干模式/定量相位衬度"的入口。
#
# 严格书约定：前向变换 ifft2*N²（+j 指数）、逆变换 fft2/N²（-j 指数）——同 ch2/ch3。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 版连通用 regionprops 定位；此处用阈值 + 质心公式（numpy 足够）。
# 2) imshow/imagesc 坐标原点差异同前几章：plt.imshow 加 origin='lower'。

import sys
import numpy as np
import matplotlib.pyplot as plt

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

# ----- 1. 参数与物体 -----
N = 512
lam = 532e-9
dx = 5e-6
L = N * dx
z_prop = 0.05
k0 = 2 * np.pi / lam
x = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(x, x)

U_obj = np.ones((N, N), dtype=complex)
radius = 40e-6                                       # 粒子半径 40 µm
mask_amp = ((X + 0.4e-3)**2 + Y**2) <= radius**2     # 左：纯振幅（强吸收、无相位）
U_obj[mask_amp] = 0.1                                 # 振幅透过率 0.1（原 MATLAB 版为 0.1，旧 .py 为 0，统一取 0.1）
mask_ph = ((X - 0.4e-3)**2 + Y**2) <= radius**2      # 右：纯相位（无吸收、延迟 1.5 rad）
U_obj[mask_ph] = np.exp(1j * 1.5)

# ----- 2. 书式傅里叶协议 + 角谱正/反向传播 -----
def ft2_book(U):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

df = 1 / L
fx = (np.arange(N) - N // 2) * df
FY, FX = np.meshgrid(fx, fx)
KX, KY = 2 * np.pi * FX, 2 * np.pi * FY
root = 1 - (KX**2 + KY**2) / k0**2
root[root < 0] = 0                                   # 滤除消逝波
H_fwd = np.exp(-1j * k0 * z_prop * np.sqrt(root))
H_bwd = np.exp(+1j * k0 * z_prop * np.sqrt(root))    # 反向传播：z 取 -z_prop

U_det = ift2_book(ft2_book(U_obj) * H_fwd)           # 传播到相机面
I_holo = np.abs(U_det)**2                            # 相机只记强度
U_rec = ift2_book(ft2_book(np.sqrt(I_holo)) * H_bwd) # sqrt 强度当复振幅估计，反向传播

# ----- 3. 定量对比：振幅通道里"看得见"谁？相位通道里"读得出"谁？ -----
amp = np.abs(U_rec)
ph = np.angle(U_rec)
bg = amp[(X**2 + Y**2) > (1.0e-3)**2].mean()         # 远场背景振幅
c_amp_obj = (bg - amp[mask_amp].mean()) / bg         # 振幅粒子在重建振幅中的对比度
c_ph_obj = abs(amp[mask_ph].mean() - bg) / bg        # 相位粒子在重建振幅中的对比度
ph_dome = np.angle(np.mean(U_rec[mask_ph]))          # 相位粒子中心读出的相位

print('--- 重建通道对比（同轴 sqrt 强度重建） ---')
print(f'振幅粒子（透过率 0.1）：重建振幅对比度 = {c_amp_obj:.3f}（清晰黑盘，形貌真实）')
print(f'相位粒子（延迟 1.5 rad，零吸收）：重建振幅对比度 = {c_ph_obj:.3f}'
      '（全部来自幽灵环干涉伪影，不是形貌）')
print(f'相位通道在相位粒子中心读出 = {ph_dome:.3f} rad（真值 1.5，被背景干涉压低、非线性）')
print('结论：振幅通道会误导（伪影像形貌）、相位通道读数失真 -> 定量相位成像'
      '必须走复全息图路线（OSH 的 ic/is 双通道，或 ch5 相干模式）')

# ----- 4. 可视化 -----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
ext = [x[0] * 1e3, x[-1] * 1e3, x[0] * 1e3, x[-1] * 1e3]
ax[0].imshow(np.abs(U_obj), extent=ext, cmap='gray', vmin=0, vmax=1.2)
ax[0].set_title('1. Object amplitude (0.1 absorber | 1.5 rad dome)')
ax[1].imshow(amp, extent=ext, cmap='gray')
ax[1].set_title('2. Reconstructed amplitude (absorber clear, phase-ghost ring)')
ax[2].imshow(ph, extent=ext, cmap='jet', vmin=-np.pi, vmax=np.pi)
ax[2].set_title('3. Reconstructed phase (1.5 rad dome, corrupted)')
for a in ax.flat:
    a.set_xlabel('x (mm)'); a.set_ylabel('y (mm)')
plt.tight_layout()
plt.savefig('sim-amplitude-phase.png', dpi=150, bbox_inches='tight')
plt.show()
