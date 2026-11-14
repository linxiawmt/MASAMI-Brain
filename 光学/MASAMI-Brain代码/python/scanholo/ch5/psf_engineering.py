# ===== 书 5.3 · PSF 工程（式 5.3-1 ~ 5.3-7）：反向曲率球面波把横向分辨率提高 2 倍 =====
# 标准 OSH（5.3a）：p1=1（球面波@w0）+ p2=delta（平面波@w0+Omega），TDFZP 扫描（式 5.3-1）。
#   针孔物体的 cos 通道全息图（式 5.3-2）：i_c ~ sin[k0 r^2 / (2 z0)]，重建像距 z0，
#   NA = r/z0，分辨率 lam0/NA。
# PSF 工程（5.3b）：把 p2 换成点源实现（点源 p 在光瞳面前 z_i 处 -> p2 = exp[-j k0 r^2/(2 z_i)]，
#   式 5.3-3），物面上得到会聚球面波（式 5.3-4：a = (jk0/2 pi z1) exp[+j k0 r^2/(2 z1)]）。
#   成像条件（式 5.3-5，原书印的就是这个形式，PDF 300 dpi 放大确认无误）：
#       1/(f+z_i) + 1/(f+z0+z1) = 1/f   <=>   f^2 = z_i (z0+z1)
#   选 z_i=2f、z0=z1=f/4 精确满足。波动光学复核：焦平面会聚波半径 R_f = f^2/z_i
#   （F{p2}|_{kx=k0 x/f} = exp[+j k0 z_i r^2/(2 f^2)]），再传 z0 后半径 R_f - z0 = z1 —— 与 (5.3-5) 等价。
#   物面上是"曲率相反的球面波对"，扫描强度（式 5.3-6）
#   = A' + B' sin[k0 r^2 / z0 + Omega t] —— 啁啾率翻倍，新全息图（式 5.3-7）
#   i_c ~ sin[k0 r^2 / z0]，重建像距 z0/2，NA = r/(z0/2)，分辨率 lam0/(2 NA)。
#
# 本脚本验证：a) 式 5.3-5 成像条件（书本数字 + 等价形式 f^2=z_i(z0+z1)）；b) 式 5.3-3 -> 5.3-4 会聚波曲率
#             R = f^2/z_i - z0；c) 两张全息图的聚焦距离；d) 焦斑 FWHM 之比 ~ 2
#             （超越孔径 Rayleigh 极限 2 倍，Indebetouw 2002 理论 / 2005 光学实验证实）。
#
# 严格书约定：前向变换 +j 指数 -> numpy 用 ifft2*N^2；逆变换 -j 指数 -> fft2/N^2。
# MATLAB / Python 差异提示：plt.imshow 加 origin='lower'（MATLAB imagesc 默认 y 向上）。

import numpy as np
import matplotlib.pyplot as plt

# ----- 0. 网格与书协议傅里叶封装 -----
# 参数设计（三个约束同时满足）：
#   a) 角谱传递函数无混叠：d < k0*N*dx^2/(2*pi) = N*dx^2/lam0（dx=2, N=1024 -> 8192）；
#   b) 菲涅尔数 N_F = r^2/(lam0*z0) ~ 4：N_F 太小时出现显著"焦移"（Li & Wolf 效应，
#      低 N_F 会聚波的轴上强度峰明显前移，见第 4 步扫描曲线的注释）；
#   c) 焦斑可采样：FWHM_std = 3.232*lam0*z0/(2*pi*r) >= ~7 px。
# 取 z0=7040, r=120: N_F=4.1, FWHM_std=15.1（7.5 px），FWHM_eng=7.5（3.8 px）。
N = 1024
dx = 2.0
lam0 = 0.5
k0 = 2 * np.pi / lam0
z0 = 7040.0                    # 物面到焦点 c 的距离
r_ap = 120.0                   # 全息图孔径半径（式 5.3-2 讨论中的 r）

x = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(x, x)
R2 = X**2 + Y**2
KX = 2 * np.pi * np.fft.fftfreq(N, d=dx)
KX2, KY2 = np.meshgrid(KX, KX)
K2 = KX2**2 + KY2**2

def fwd(u):                                   # 书前向变换（+j 指数）
    return np.fft.ifft2(u) * N**2

def inv(U):                                   # 书逆变换（-j 指数）
    return np.fft.fft2(U) / N**2

def propagate(u, d):                          # 书式 (2.3-13) 传递函数，向 +z 传播 d
    return inv(fwd(u) * np.exp(-1j * k0 * d) * np.exp(1j * d * K2 / (2 * k0)))

ap = R2 <= r_ap**2             # 全息图孔径（半径 r）
row = N // 2

# ----- 1. 式 5.3-5 成像条件（原书即印 1/(f+z_i)+1/(f+z0+z1)=1/f，无误）-----
# 交叉验证三条路：a) 书本数字代入；b) 等价形式 f^2 = z_i(z0+z1)（交叉相乘可证）；
# c) 波动光学链：R_f = f^2/z_i（焦平面）-> R = R_f - z0 = z1（物面），与 b) 恒等。
# （若把 p2 的下标误读成 z1，则条件变成 f^2 = z1(z0+z1)，书数 z_i=2f, z0=z1=f/4 不满足——
#   故 p2 下标必为 z_i：它由"点源放在光瞳面前 z_i 处"实现，见书 p.146 的实现说明。）
f = 8000.0                     # 透镜焦距
zi, z0c, z1c = 2 * f, f / 4, f / 4   # 书中选值：z_i=2f, z0=z1=f/4
lhs = 1 / (f + zi) + 1 / (f + z0c + z1c)          # 式 (5.3-5) 左端
print(f'--- 式 5.3-5 成像条件核验（f={f:.0f}, z_i=2f, z0=z1=f/4）---')
print(f'1/(f+z_i)+1/(f+z0+z1) = {lhs * f:.6f} /f (应 = 1.000000，书本数字精确成立)')
print(f'等价形式：z_i = f^2/(z0+z1) = {f**2 / (z0c + z1c) / f:.4f} f (应 = 2f，一致)')

# ----- 2. 核验式 5.3-3 -> 5.3-4：物面会聚波的曲率半径 = f^2/z_i - z0（= z1）-----
# 点源（p2 的实现）经透镜成像：焦平面场是会聚 chirp，半径 f^2/z_i（式 5.3-3 的
# F{p2}|_{kx=k0 x/f}）；再卷积 h(x,y;z0) 传到物面（式 5.3-3），应得会聚波半径 z1。
R_f = f**2 / zi                                  # 焦平面会聚波半径（本参数下 = 4000）
u_f = np.exp(1j * k0 * R2 / (2 * R_f))           # 焦平面场（会聚）
a_field = propagate(u_f, z0c)                    # 式 (5.3-3)：再卷积 h(x,y;z0)

seg = slice(N // 2 - 60, N // 2 + 60)            # 只在采样充分的中心区拟合
ph = np.unwrap(np.angle(a_field[row, seg]))
rho = x[seg]
cf = np.polyfit(rho**2, ph, 1)                   # 相位 ~ cf[0] * r^2 + cf[1]
R_num = k0 / (2 * cf[0])
print('--- 式 5.3-4 会聚波曲率核验（点源实现 z_i=2f）---')
print(f'理论 R = f^2/z_i - z0 = {f**2 / zi - z0c:.1f}  (即 z1 = {z1c:.1f})')
print(f'数值拟合 R_num        = {R_num:.1f}  (相对误差 {abs(R_num - (f**2 / zi - z0c)) / (f**2 / zi - z0c) * 100:.2f}%)')

# ----- 3. 两张针孔全息图（式 5.3-2 / 5.3-7 的复全息图版，取 +chirp 支）-----
H_std = ap * np.exp(1j * k0 * R2 / (2 * z0))     # 标准 TDFZP 记录，重建像距 z0
H_eng = ap * np.exp(1j * k0 * R2 / z0)           # 反向曲率对记录，重建像距 z0/2
s_std = ap * np.sin(k0 * R2 / (2 * z0))          # 式 (5.3-2) 的 cos 通道（显示用）
s_eng = ap * np.sin(k0 * R2 / z0)                # 式 (5.3-7)

d_grid = np.linspace(0.4 * z0, 1.1 * z0, 29)     # 聚焦距离扫描（覆盖 z0/2 与 z0，且 < AS 采样上限）
sharp_std = np.array([np.abs(propagate(H_std, d)).max()**2 for d in d_grid])
sharp_eng = np.array([np.abs(propagate(H_eng, d)).max()**2 for d in d_grid])
d_best_std = d_grid[np.argmax(sharp_std)]
d_best_eng = d_grid[np.argmax(sharp_eng)]
print('--- 聚焦距离扫描 ---')
print(f'标准全息图最锐距离 = {d_best_std:.0f} (几何焦点 z0 = {z0:.0f}；N_F={r_ap**2 / (lam0 * z0):.1f} '
      '的低菲涅尔数焦移使峰值略前移，Li & Wolf 效应)')
print(f'工程全息图最锐距离 = {d_best_eng:.0f} (几何焦点 z0/2 = {z0 / 2:.0f})')

# ----- 4. 焦斑 FWHM（线性插值找半高交点）-----
def fwhm(profile, step):
    I = profile / profile.max()
    half = 0.5
    pk = int(np.argmax(I))
    hi = pk
    while hi < len(I) - 1 and I[hi + 1] >= half:
        hi += 1
    xr = hi + (I[hi] - half) / (I[hi] - I[hi + 1]) if hi + 1 < len(I) else hi
    lo = pk
    while lo > 0 and I[lo - 1] >= half:
        lo -= 1
    xl = lo - (I[lo] - half) / (I[lo] - I[lo - 1]) if lo - 1 >= 0 else lo
    return (xr - xl) * step

U_std = propagate(H_std, z0)       # 标准焦斑（z0 处）
U_eng = propagate(H_eng, z0 / 2)   # 工程焦斑（z0/2 处）
w_std, w_eng = fwhm(np.abs(U_std[row, :])**2, dx), fwhm(np.abs(U_eng[row, :])**2, dx)
print('--- 焦斑 FWHM（横向分辨率）---')
print(f'标准：FWHM = {w_std:.2f}（Airy 理论 {3.232 * z0 / (k0 * r_ap):.2f}）—— 分辨率 lam0/NA，NA = r/z0')
print(f'工程：FWHM = {w_eng:.2f}（Airy 理论 {3.232 * (z0 / 2) / (k0 * r_ap):.2f}）—— 分辨率 lam0/(2NA)')
print(f'FWHM 之比 = {w_std / w_eng:.3f}（理论 2.0：超越孔径 Rayleigh 极限 2 倍）')

# ----- 5. 可视化 -----
fig, ax = plt.subplots(2, 2, figsize=(12.5, 10))
ext = [x[0], x[-1], x[0], x[-1]]
ax[0, 0].imshow(s_std, extent=ext, cmap='gray')
ax[0, 0].set_title('1. Std hologram $i_c\\sim\\sin(k_0r^2/2z_0)$  (Eq.5.3-2), r=%.0f' % r_ap)
ax[0, 0].set_xlim(-250, 250); ax[0, 0].set_ylim(-250, 250)
ax[0, 1].imshow(s_eng, extent=ext, cmap='gray')
ax[0, 1].set_title('2. Engineered $i_c\\sim\\sin(k_0r^2/z_0)$  (Eq.5.3-7, chirp x2)')
ax[0, 1].set_xlim(-250, 250); ax[0, 1].set_ylim(-250, 250)
I_std = np.abs(U_std[row, :])**2; I_eng = np.abs(U_eng[row, :])**2
ax[1, 0].plot(x, I_std / I_std.max(), label='std @z0, FWHM=%.1f' % w_std)
ax[1, 0].plot(x, I_eng / I_eng.max(), label='eng @z0/2, FWHM=%.1f' % w_eng)
ax[1, 0].set_xlim(-30, 30); ax[1, 0].legend(fontsize=9); ax[1, 0].grid(alpha=0.3)
ax[1, 0].set_title('3. Focal spot cross-sections (FWHM ratio ~ 2)')
ax[1, 0].set_xlabel('x'); ax[1, 0].set_ylabel('normalized intensity')
ax[1, 1].plot(d_grid / z0, sharp_std / sharp_std.max(), 'o-', ms=3, label='std (peak @ z0)')
ax[1, 1].plot(d_grid / z0, sharp_eng / sharp_eng.max(), 's-', ms=3, label='eng (peak @ z0/2)')
ax[1, 1].legend(fontsize=9); ax[1, 1].grid(alpha=0.3)
ax[1, 1].set_title('4. Focus sharpness vs reconstruction distance')
ax[1, 1].set_xlabel('reconstruction distance d / z0'); ax[1, 1].set_ylabel('peak |U|^2 (norm.)')
plt.tight_layout()
plt.savefig('sim-psf-engineering.png', dpi=150, bbox_inches='tight')
plt.show()
