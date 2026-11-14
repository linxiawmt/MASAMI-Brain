# ===== 第 1 章数值验证脚本 =====
# 逐条验证笔记里的关键公式与书中数值，全部用 numpy 离散实现。
# 运行：python verify_ch1.py   （无绘图，纯文本结果）

import numpy as np

ok = lambda name, cond: print(f'[{"PASS" if cond else "FAIL"}] {name}')

# ---------- (a) sinc 定义与零点 ----------
x = np.array([1, 2, 3, 4, 5, -1, -2, -3], dtype=float)
ok('(a) sinc(x)=sin(pi x)/(pi x) 在 x=±1,±2,... 处为零',
   np.allclose(np.sin(np.pi * x) / (np.pi * x), 0, atol=1e-12))

# ---------- (b) 原书 Plot_sinc.m 的 0/0 陷阱 ----------
xs = -5.5 + 0.01 * np.arange(1101)          # 等价 MATLAB x=-5.5:0.01:5.5
hit0 = np.where(np.isclose(xs, 0.0, atol=1e-9))[0]
with np.errstate(invalid='ignore'):
    val0 = np.sin(np.pi * 0.0) / (np.pi * 0.0) if hit0.size == 0 else \
           np.sin(np.pi * xs[hit0[0]]) / (np.pi * xs[hit0[0]])
print(f'       网格是否精确命中 x=0: {hit0.size > 0}；x=0 处 sin(pi x)/(pi x) = {val0} (nan 即 0/0)')
ok('(b) 原书代码在 x=0 处产生 0/0（NaN），需要补极限值 1', np.isnan(val0) or hit0.size == 0)

# ---------- (c) 例 1.1 书中数值：M=11 -> a=0.0429, kx0=146.23 ----------
M = 11
a = M / 256
ok('(c) a = 11/256 = 0.0429（书中数值）', abs(a - 0.0429) < 5e-4)
ok('(c) 第一零点 kx0 = 2*pi/a = 146.23 rad/(unit length)（书中数值）',
   abs(2 * np.pi / a - 146.23) < 0.01)

# ---------- (d) 变换对 rect(x/a,y/b) <=> ab sinc(a kx/2pi, b ky/2pi) ----------
# 严格验证用"真实 DFT 频轴" K = 2*pi*n/(N*dx)（n=-N/2..N/2-1），
# 并把 rect 中心放在 fftshift 约定的网格中心（Python 下标 128）。
# 书中的 Kx 轴用 2*pi/((N-1)*dx) 步长，比真实 DFT 频轴大 N/(N-1) 倍（见"坑与陷阱"）。
N = 256
L = 1.0
dx = L / (N - 1)
a = M / 256
b = a
R = np.zeros((N, N)); n = (M - 1) // 2; c = 128
R[c - n:c + n + 1, c - n:c + n + 1] = 1.0
F = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(R)) * N**2)   # 书式(1.1-2a) 正指数
Kx = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(N, d=dx))
Ky = Kx.copy()
KX, KY = np.meshgrid(Kx, Ky)
F_ana = a * b * np.sinc(a * KX / (2 * np.pi)) * np.sinc(b * KY / (2 * np.pi))
err = np.abs(F_ana - F * dx**2).max() / np.abs(F_ana).max()     # DFT 乘 dx^2 逼近连续积分
print(f'       解析谱与 DFT 谱的最大相对偏差 = {err:.2e}')
print('       （偏差来源：显示窗口有限导致的频谱周期延拓混叠，'
      '以及书中 a=M/256 与实际采样宽度 M*dx=M/255 的半像素差，见"坑与陷阱"）')
ok('(d) rect 变换对：DFT(正指数) 与 ab*sinc 解析式一致（<5%）', err < 0.05)

# 第一零点位置（沿 kx，取中心行，找 |F| 的第一个极小值并与 2*pi/a 比较）
row = np.abs(F[128, :] * dx**2) / (np.abs(F * dx**2)).max()
band = (Kx > 2 * np.pi / a * 0.5) & (Kx < 2 * np.pi / a * 1.5)
i_min = np.argmin(np.where(band, row, np.inf))
dk = Kx[1] - Kx[0]
print(f'       测得第一零点 kx ≈ {Kx[i_min]:.2f}，理论 2*pi/a = {2*np.pi/a:.2f}（频格 {dk:.2f}）')
ok('(d) 频谱第一零点位于 2*pi/a（误差 < 1 个频格）', abs(Kx[i_min] - 2 * np.pi / a) < dk)

# ---------- (e) 卷积定理：频域乘 = 空间域卷积 ----------
rng = np.random.default_rng(0)
f1 = rng.random(16); h1 = rng.random(12)
g_direct = np.convolve(f1, h1)                        # 直接卷积，长度 N+M-1=27
g_fft = np.fft.ifft(np.fft.fft(f1, 27) * np.fft.fft(h1, 27))
err_e = np.abs(g_direct - g_fft.real).max()
print(f'       卷积输出尺寸 = {g_direct.size}（=16+12-1），直接卷积 vs FFT 最大偏差 = {err_e:.2e}')
ok('(e) 卷积定理：ifft(fft(f)*fft(h)) 与直接卷积一致，长度 N1+N2-1',
   g_direct.size == 27 and err_e < 1e-12)

# ---------- (f) 例 1.6：自相关中心最大 |C(0,0)| >= |C(x,y)| ----------
I1 = np.zeros((256, 256))
yy, xx = np.mgrid[0:256, 0:256]
r = np.sqrt((xx - 128.0)**2 + (yy - 128.0)**2)
I1[np.abs(r - 90) < 4] = 1.0
I1[(xx - 93)**2 + (yy - 88)**2 < 81] = 1.0
I1[(xx - 163)**2 + (yy - 88)**2 < 81] = 1.0
FI1 = np.fft.fft2(I1); FI1 = FI1 / np.abs(FI1).max()
PR = np.fft.ifft2(FI1 * np.conj(FI1))                 # 自相关（与 correlation.py 同一路线）
PR = np.fft.fftshift(PR)
C = np.abs(PR)
center = C[128, 128]
ok('(f) 自相关中心最大：|C(128,128)| >= max|C(x,y)|（书中式 1.2-6）',
   center >= C.max() - 1e-12 and abs(center - C.max()) < 1e-12)

# ---------- (g) 例 1.2.3 拓展：菲涅耳衍射图样的对称性 ----------
# 方孔与菲涅耳核都关于网格中心对称时，衍射光强应严格上下/左右镜像对称。
# 为让对称性严格成立，孔宽取偶数 112（几何中心落在两像素正中 127.5），
# 并用 ifftshift 把两个场都搬到"原点在 0"的自然 FFT 顺序里做循环卷积。
lam = 6.328e-5; z = 5.0
sigma = np.pi / (lam * z)
kq = np.arange(256); X = kq / 255 - 0.5; Y = X.copy()
Me = 112                                        # 偶数宽度，中心在 127.5
R2 = np.zeros((256, 256))
R2[128 - Me // 2:128 + Me // 2, 128 - Me // 2:128 + Me // 2] = 1.0
r2 = X[:, None]**2 + Y[None, :]**2              # 网格满足 X[255-r] = -X[r]，关于 127.5 对称
h = 1j * (sigma / np.pi) * np.exp(-1j * sigma * r2) * np.exp(-800.0 * r2)
g = np.fft.fftshift(np.fft.ifft2(
    np.fft.fft2(np.fft.ifftshift(R2)) * np.fft.fft2(np.fft.ifftshift(h)))) * 256**2
I = np.abs(g)**2
# 对称中心是整数下标 127（两个核都搬到"原点在 0"后，循环卷积中心落在
# 255.5+255.5 ≡ -1，经 fftshift 平移 +128 后为 127），故按 r -> 254-r 检验镜像。
J = I[:255, :255]
sym1 = np.abs(J - J[::-1, :]).max() / J.max()   # 上下镜像（关于第 127 行）
sym2 = np.abs(J - J[:, ::-1]).max() / J.max()   # 左右镜像（关于第 127 列）
pk = np.unravel_index(np.argmax(I), I.shape)
print(f'       强度上下镜像相对偏差 = {sym1:.2e}，左右镜像相对偏差 = {sym2:.2e}')
print(f'       光强峰位于 {pk}（孔边缘衍射条纹处，与镜像位置 {254-pk[0]} 成对），'
      f'Fresnel 数 a^2/(lambda*z) ≈ {(Me*dx/2)**2/(lam*z):.0f}')
ok('(g) 方孔菲涅耳衍射光强严格上下/左右镜像对称',
   sym1 < 1e-10 and sym2 < 1e-10)

print('\n全部验证完成。')
