# ===== 书 4.2 节 · 3-D 全息电视的三大瓶颈（式 4.2-6 ~ 4.2-10 + 数据率） =====
# A. 空间频率分辨率 -> 视角：
#    点源全息图 sin[(k0/2z0)(x²+y²)] 的局部空间频率 (4.2-6)：
#        f_local = (1/2π) d/dx[(k0/2z0)x²] = x/(λ0·z0)
#    全息图孔径半宽 x_max 处最大 (4.2-7)：f_max = x_max/(λ0·z0)
#    几何 (4.2-8)：NA = sin(θ/2) = x_max/z0  ->  (4.2-9)：NA = λ0·f_max = λ0·f0
#    即：SLM 分辨率 f0 直接决定视角 θ = 2·arcsin(λ0·f0)。
# B. 空间分辨率（像素数）：Nyquist 采样 fs = 2·f0；l×l 全息图需
#        N = (l·fs)² = (2l·f0)² = (2l·NA/λ0)²                (4.2-10)
# C. 数据传输：N×8bit×30帧/s；HPO（只保留水平视差）可把像素数砍到 512×(2l·NA/λ0)。
#
# 本脚本复核书中的全部数值，并验证书 4.2-C 的一处笔误：
#    书印 HPO 像素公式 512×(2l×2NA/λ0) —— 代入得 32.4 million，与书自己给的
#    "approximately 17 million" 矛盾；正确式应为 512×(2l·NA/λ0) ≈ 16.2 million（见笔记勘误）。
#
# MATLAB / Python 差异提示：纯数值脚本，无傅里叶约定问题。

import sys
import numpy as np
import matplotlib.pyplot as plt

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

lam0 = 0.6328e-3          # mm（He-Ne）
l = 20.0                  # 全息图边长 mm（20mm × 20mm）
theta_full = 60.0         # 视角（度）
NA = np.sin(np.radians(theta_full / 2))

# ----- A. 视角 vs SLM 分辨率（Table 4.1） -----
f0_list = np.array([8.0, 100.0, 500.0])           # lp/mm
theta_list = 2 * np.degrees(np.arcsin(lam0 * f0_list))
print('--- A. 视角（Table 4.1，λ0=0.6328 µm） ---')
for f0, th, dev in zip(f0_list, theta_list, ['EBSLM/Hamamatsu', 'PALSLM/Hamamatsu', 'not available']):
    print(f'    f0 = {f0:5.0f} lp/mm -> θ = {th:5.2f} deg   （书 Table 4.1：'
          f'{["0.6", "6.8", "34"][list(f0_list).index(f0)]}°，{dev}）')
print('    离轴全息 45° 记录角 -> 载频 ~1000 lp/mm：远超现有 SLM（同轴 OSH 才上得了 SLM）')

# ----- B. 像素数（4.2-10） -----
f0_req = NA / lam0                                    # 满足 60° 视角所需的 SLM 分辨率（4.2-9 反推）
fs = 2 * f0_req                                       # Nyquist：fs = 2 f0
N_pix = (l * fs)**2
print('--- B. 像素数（式 4.2-10，l=20mm，θ=60°） ---')
print(f'    所需 f0 = NA/λ0 = {f0_req:.0f} lp/mm，fs = 2f0 = {fs:.0f} lp/mm')
print(f'    N = (l·fs)² = (2l·NA/λ0)² = {N_pix/1e9:.2f} billion（书：about 1.1 billion，'
      '按 λ0=0.6 µm 计恰为 1.1）')
print('    对照 Canon D60：3072×2048 ≈ 6.3 MP（67.7 lp/mm，7.4 µm 像素）——差 160 倍')

# ----- C. 数据率与 HPO -----
N_book = 1.1e9                                        # 书采用 1.1 billion
rate_full = N_pix * 8 * 30 / 1e12
print('--- C. 数据率 ---')
print(f'    全视差：{N_pix:.3e} samples × 8 bit × 30 fps = {rate_full:.2f} Tbit/s'
      f'（书用 1.1 billion -> 0.26 Tbit/s）')
n_hpo_book = 512 * (2 * l * 2 * NA / lam0)            # 书印公式（多一个因子 2）
n_hpo_fix = 512 * (2 * l * NA / lam0)                 # 与 17 million 自洽的式
print(f'    HPO（512 条垂直线）：')
print(f'      书印式 512×(2l×2NA/λ0) = {n_hpo_book/1e6:.1f} million  <- 与书自己的 "17 million" 矛盾（笔误）')
print(f'      修正式 512×(2l·NA/λ0)  = {n_hpo_fix/1e6:.1f} million ≈ 17 million ✓'
      f'  -> {n_hpo_fix*8*30/1e9:.1f} Gbit/s（书：4 Gbit/s）')
print('    光纤 40 Gbit/s 可达 -> 实时全息电视可行，前提是 HPO 电子全息记录（HPO-OSH，1-D TDFZP 扫描）')

# ----- 可视化（对应笔记图 sim-slm-budget.png） -----
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
f0s = np.linspace(1, 600, 400)
th = 2 * np.degrees(np.arcsin(np.clip(lam0 * f0s, 0, 1)))
ax[0].plot(f0s, th, 'b-', lw=1.5)
for f0, t0, name in [(8, 0.58, 'EBSLM'), (100, 7.27, 'PALSLM'), (500, 36.9, 'ideal')]:
    ax[0].plot(f0, 2 * np.degrees(np.arcsin(lam0 * f0)), 'ro')
    ax[0].annotate(name, xy=(f0, 2 * np.degrees(np.arcsin(lam0 * f0))),
                   xytext=(f0 * 0.55, 2 * np.degrees(np.arcsin(lam0 * f0)) + 6), fontsize=9)
ax[0].axhline(60, color='r', ls='--', lw=0.8, label='60° target')
ax[0].set_xlabel("SLM resolution $f_0$ (lp/mm)"); ax[0].set_ylabel('viewing angle θ (deg)')
ax[0].set_title('Viewing angle vs SLM resolution (Eq. 4.2-9)')
ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)

rates = [N_book * 8 * 30 / 1e9, n_hpo_fix * 8 * 30 / 1e9]
ax[1].bar(['full parallax', 'HPO (512 lines)'], rates, color=['#c0392b', '#27ae60'], width=0.5)
ax[1].set_yscale('log'); ax[1].set_ylabel('serial data rate (Gbit/s)')
ax[1].set_title('Data rate: 8 bit × 30 fps')
for i, rv in enumerate(rates):
    ax[1].text(i, rv * 1.2, f'{rv:.0f} Gbit/s', ha='center', fontsize=9)
ax[1].axhline(40, color='k', ls=':', label='fiber 40 Gbit/s')
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('sim-slm-budget.png', dpi=150, bbox_inches='tight')
plt.show()
