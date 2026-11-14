# ===== 书 Example 4.1 · 全息放大率（式 4.2-1 ~ 4.2-5） =====
# 三点物体 (4.2-1)：δ(x,y,z-z0) + δ(x-x0,y,z-z0) + δ(x,y,z-(z0+Δz0))
# 正弦全息图 (4.2-2)：
#   H_3p(x,y) ~ sin[(k0/2z0)(x²+y²)] + sin{(k0/2z0)[(x-x0)²+y²]} + sin[(k0/2(z0+Δz0))(x²+y²)]
# a) 全息图缩放 M（(4.2-3)/(4.2-4)）：H(Mx,My) 中各项变为距离 z/M² 的 FZP、
#    横向位置 x0/M -> M_lat = 1/M, M_long = 1/M²。
# b) 波长缩放 m（(4.2-5)）：重建波长 λ0' = m·λ0（波数 k0/m）-> 焦距 z -> z0/m（m>1 压缩），
#    M_lat = 1, M_long = 1/m。
# c) 无畸变条件：m = 1/M -> M_lat = M_long（Gabor 1949 的原始想法）。
#
# 数值验证（本笔记实跑，结论）：
#   a) M = 0.5（放大 2 倍）：缩放全息图在 z0/M² = 120 mm 聚焦出两点，横向间距 = x0/M = 2.00 mm；
#      轴上深度响应两个峰 120/200 mm -> 纵向间距 Δz0/M² = 80 mm = 4·Δz0（M_long = 1/M²）✓
#   b) m = 2（重建波数 k0/m）：原始全息图在 z0/m = 50 mm 聚焦，横向位置不变（M_lat = 1）✓
#
# 采样与混叠账（ch3 坑的延续，两组实验用不同 z0 的原因）：
#   FZP 条纹最细周期（视场边缘）= λ0·z0/x_max，需 > 2·dx：
#     a) 缩放全息图等效周期 = λ0·z0/(M²·x_max) -> z0 = 30 mm 即满足；
#     b) 原始全息图 -> 需 z0 >= 2·dx·x_max/λ0 ≈ 81 mm，故取 100 mm。
#   重建传递函数 e^{+jz'k²/2k'} 的相邻格点相位跳变 z'·k·Δk/k' < π：
#     深度扫描全程固定窗（孔径不随 z' 变化，否则深度响应被窗趋势污染），
#     窗缘 12 次方软滚降（硬圆孔引入衍射纹波）。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 1-based、Python 0-based：中心像素 MATLAB (N/2+1) <-> Python N//2。
# 2) plt.imshow 加 origin='lower' 对齐数学坐标。

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, signal

if hasattr(sys.stdout, 'reconfigure'):      # Windows 控制台中文输出
    sys.stdout.reconfigure(encoding='utf-8')

# ----- 1. 公共参数 -----
lam = 0.6328e-6                 # He-Ne 波长 (m)
k0 = 2 * np.pi / lam
x0 = 1.0e-3                     # 前两点横向间距 (m)
N = 512
dx = 10e-6                      # 全息图像素 (m)
M = 0.5                         # a) 全息图缩放因子（<1 即放大 1/M 倍）
m = 2.0                         # b) 波长缩放因子（重建波数 k0/m）

x = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(x, x)
x_max = x[-1]

def ft2_book(U):
    # 书式前向变换（+j 指数）
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(U))) * N**2

def ift2_book(A):
    # 书式逆变换（-j 指数）
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(A))) / N**2

dk = 2 * np.pi / (N * dx)             # 角频率格距 (rad/m)
kx1d = (np.arange(N) - N // 2) * dk
KY, KX = np.meshgrid(kx1d, kx1d)
KR = np.hypot(KX, KY)
KR2 = KX**2 + KY**2
k_nyq = 0.55 * np.pi / dx             # 通用窗上限

def recon_transfer(zp, k=k0, k_win=None):
    # 与 h*(z') = +j(k/2πz')e^{-j k r²/2z'} 对应的频域传递（= conj(OTFosh)，见 ch3 式 3.5-2a）
    # 圆窗压二次相位的 Poisson 混叠：相邻格点跳变 z'·k·Δk/k' < π（留 20% 裕量），
    # 窗缘 12 次方软滚降。
    if k_win is None:
        k_win = min(k_nyq, 0.8 * np.pi * k / (zp * dk))
    return (1.0 / (1.0 + (KR / k_win)**12)) * np.exp(1j * zp * KR2 / (2 * k))

def recon(Hspec, zp, k=k0, k_win=None):
    rec = np.abs(ift2_book(Hspec * recon_transfer(zp, k, k_win)))**2
    return rec / rec.max()

def peaks_of(rec, thr=0.5):
    mask = (ndimage.maximum_filter(rec, size=31) == rec) & (rec > thr)
    lab, nlab = ndimage.label(mask)
    pts = sorted(((p[1] - N / 2) * dx * 1e3, (p[0] - N / 2) * dx * 1e3, rec[p[0], p[1]])
                 for p in ndimage.maximum_position(rec, lab, range(1, nlab + 1)))
    return pts

def report(rec, label, thr=0.5):
    pts = peaks_of(rec, thr)
    print(f'{label}：{len(pts)} 个峰（阈值 {thr}），'
          + ', '.join(f'({px:+.2f},{py:+.2f},I={ii:.2f})' for px, py, ii in pts[:4]))
    return pts

def make_h3p(z0, dz0):
    # 三点物体的正弦全息图 (4.2-2)（不加常数偏置）
    return (np.sin(k0 * (X**2 + Y**2) / (2 * z0))
            + np.sin(k0 * ((X - x0)**2 + Y**2) / (2 * z0))
            + np.sin(k0 * (X**2 + Y**2) / (2 * (z0 + dz0))))

# ================= a) 全息图缩放 M = 0.5 =================
z0a, dz0a = 0.030, 0.020        # 物体深度（缩放全息图等效条纹周期 λ0z0/(M²x_max) 采样允许更近）
Hs = (np.sin(k0 * M**2 * (X**2 + Y**2) / (2 * z0a))                       # (4.2-4) 逐项解析代入
      + np.sin(k0 * M**2 * ((X - x0 / M)**2 + Y**2) / (2 * z0a))
      + np.sin(k0 * M**2 * (X**2 + Y**2) / (2 * (z0a + dz0a))))
FHs = ft2_book(Hs)

print('--- a) 全息图缩放 M = 0.5（理论：M_lat = 1/M = 2, M_long = 1/M^2 = 4） ---')
zs_scan = np.linspace(0.08, 0.26, 37)
k_win_fix = min(k_nyq, 0.8 * np.pi * k0 / (zs_scan[-1] * dk))   # 深度扫描全程同一窗
r1 = recon(FHs, z0a / M**2, k_win=k_win_fix)
pts1 = report(r1, f"缩放全息图 @ z0/M² = {z0a/M**2*1e3:.0f} mm")
print(f'    两点横向间距 = {pts1[-1][0]-pts1[0][0]:.2f} mm（理论 x0/M = {x0/M*1e3:.2f} mm）')

xs_on = N // 2
I_on = np.array([recon(FHs, zp, k_win=k_win_fix)[N // 2, xs_on] for zp in zs_scan])   # 焦点在轴上中心像素，取该像素（列最大会被离焦环纹抬高）
above = I_on > 0.95 * I_on.max()          # 聚焦平台（宽度 = 该孔径的景深 ~ ±20 mm）
edges = np.diff(above.astype(int))
starts = list(zs_scan[1:][edges == 1]); ends = list(zs_scan[1:][edges == -1])
if above[0]: starts = [zs_scan[0]] + starts
if above[-1]: ends = ends + [zs_scan[-1]]
platau = [(s, e) for s, e in zip(starts, ends) if e - s >= 3 * (zs_scan[1] - zs_scan[0])]   # 滤掉纹波碎片
centers = [0.5 * (s + e) for s, e in platau]
print('    轴上深度响应聚焦平台（I>0.95·max）：'
      + ', '.join(f'{s*1e3:.0f}~{e*1e3:.0f} mm' for s, e in platau)
      + f'（理论焦深 z0/M² = {z0a/M**2*1e3:.0f}，(z0+Δz0)/M² = {(z0a+dz0a)/M**2*1e3:.0f}）')
if len(centers) >= 2:
    print(f'    平台中心间距 = {(centers[-1]-centers[0])*1e3:.0f} mm'
          f'（理论 Δz0/M² = {dz0a/M**2*1e3:.0f} mm = {1/M**2:.0f}·Δz0；偏差来自 ±20 mm 孔径景深）')

# ================= b) 波长缩放 m = 2 =================
z0b, dz0b = 0.100, 0.030        # 原始全息图边缘条纹采样要求 z0 >= 81 mm
H3p = make_h3p(z0b, dz0b)
FH3p = ft2_book(H3p)

print('--- b) 波长缩放 m = 2（重建波数 k0/m，理论：M_lat = 1, 焦深 z -> z0/m） ---')
r2 = recon(FH3p, z0b / m, k=k0 / m)
pts2 = report(r2, f'原始全息图 @ z0/m = {z0b/m*1e3:.0f} mm (m=2, 波长加倍)')
r3 = recon(FH3p, z0b)
pts3 = report(r3, f'原始全息图 @ z0 = {z0b*1e3:.0f} mm (m=1, 参考)')
print(f'    横向位置不变：{[f"{px:+.2f}" for px, py, _ in pts2[:2]]} mm（m=2）'
      f' vs {[f"{px:+.2f}" for px, py, _ in pts3[:2]]} mm（m=1）-> M_lat = 1')
print(f'    第三点焦深：(z0+Δz0)/m = {(z0b+dz0b)/m*1e3:.0f} mm (m=2) vs (z0+Δz0)/m'
      f' = {(z0b+dz0b)*1e3:.0f} mm (m=1)，纵向压缩 1/m（M_long = 1/m）')

# ================= c) 无畸变条件核对 =================
print('--- c) 无畸变条件（Gabor）：m = 1/M -> M_lat = M_long ---')
Mcheck = 1.0 / m
print(f'    M = 1/m = {Mcheck}: M_lat = 1/M = {1/Mcheck:.2f}, '
      f'M_long = 1/(M²m) = {1/(Mcheck**2*m):.2f} -> 相等 ✓')

# ----- 可视化（对应笔记图 sim-magnification.png） -----
fig, ax = plt.subplots(2, 3, figsize=(13.5, 7.5))
ext = [x[0] * 1e3, x[-1] * 1e3, x[0] * 1e3, x[-1] * 1e3]
ax[0, 0].imshow(Hs, cmap='gray', extent=ext, vmax=2.2, vmin=-1.2)
ax[0, 0].set_title(f'scaled hologram H(Mx,My), M={M} (4.2-3)')
ax[0, 1].imshow(r1**0.4, extent=ext, cmap='gray')
ax[0, 1].set_title(f"recon @ z0/M²={z0a/M**2*1e3:.0f} mm: 2 spots, dx=x0/M={x0/M*1e3:.1f} mm")
ax[0, 2].plot(zs_scan * 1e3, I_on / I_on.max(), 'b-', label='on-axis intensity')
for s, e in platau:
    ax[0, 2].axvspan(s * 1e3, e * 1e3, color='b', alpha=0.12)
ax[0, 2].axvline(z0a / M**2 * 1e3, color='r', ls='--', lw=0.8)
ax[0, 2].axvline((z0a + dz0a) / M**2 * 1e3, color='r', ls='--', lw=0.8)
ax[0, 2].set_xlabel("reconstruction depth z' (mm)"); ax[0, 2].set_ylabel('normalized intensity')
ax[0, 2].set_title('depth response of scaled hologram (M_long=1/M²)')
ax[0, 2].legend(fontsize=8); ax[0, 2].grid(alpha=0.3)
ax[1, 0].imshow(H3p, cmap='gray', extent=ext, vmax=2.2, vmin=-1.2)
ax[1, 0].set_title('H_3p(x,y): sine hologram (4.2-2), z0=100 mm')
ax[1, 1].imshow(r2**0.4, extent=ext, cmap='gray')
ax[1, 1].set_title(f"recon @ z0/m={z0b/m*1e3:.0f} mm with k'=k0/m (M_lat=1)")
ax[1, 2].imshow(r3**0.4, extent=ext, cmap='gray')
ax[1, 2].set_title(f'recon @ z0={z0b*1e3:.0f} mm with k0 (reference)')
for a in ax.flat:
    a.set_xlabel('x (mm)'); a.set_ylabel('y (mm)')
plt.tight_layout()
plt.savefig('sim-magnification.png', dpi=150, bbox_inches='tight')
plt.show()
