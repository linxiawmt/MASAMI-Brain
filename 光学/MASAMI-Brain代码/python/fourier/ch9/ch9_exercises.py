# -*- coding: utf-8 -*-
"""书 9.4 节习题 9.1-9.6 的 Python 数值实现(与 matlab/fourier/ch9/ch9_exercises.m 对应)。

  9.1 矩形线型:|gamma| = |sinc(dnu*tau)| = 0.8270,改 pc_temp 线型后与解析剖面对比
  9.2 条纹位移:相干双束、deld = 0, lam0/4, lam0/2, 3lam0/4 的条纹平移与测量局限
  9.3 可见度:A1=0.75、A2=0.20 的 V 及"提议"的解析式 vs 仿真
  9.4 空间相干 l_cr = 15/8/4/1 mm:|mu_12|、sigma_f > L/2 的失效模式与补救
  9.5 单个 5mm x 2mm 矩形孔:l_cr 从 20 mm 到 1 mm 的夫琅禾费图样演变
  9.6 高斯谢尔模型(GSM)光束:随机相位屏 + propTF 传播 1 km vs 式 (9.32)/(9.33)

运行产物:out_fig9_ex91..ex96.png + 控制台数值验证。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "ch5"))
from base import coords, circle, rect                    # noqa: E402
from ch9_temporal import pc_temp, analytic_two_beam, gamma_temporal, \
    fringe_visibility                                    # noqa: E402
from ch9_spatial import pc_spatial, mu12_exact920        # noqa: E402
from ch5_prop import propTF                              # noqa: E402  书 5.1 节传播器

OUT = []

# ============================================================
# 习题 9.1:矩形线型(式 9.29/9.30)
# ============================================================
print("=" * 62)
print("习题 9.1 矩形线型")
lam0, c, delnu, deld = 650e-9, 3e8, 2e9, 5e-2
tau = deld / c
tau_c = 1 / delnu                                    # 式 (9.30)
g_abs = abs(np.sinc(delnu * tau))                    # |gamma| = |sinc(dnu*tau)|
print(f"tau_c = 1/dnu = {tau_c:.2e} s(高斯线型是 0.664/dnu)")
print(f"|gamma(deld/c)| = |sinc(dnu*tau)| = {g_abs:.4f}"
      f"(dnu*tau = {delnu*tau:.4f};np.sinc 是归一化定义 sin(pi x)/(pi x))")
I_r, x2, (nu_arr, S_arr, dnu, lam0_, f_, w_, dels_) = pc_temp(lineshape="rect")
prof_r = I_r[I_r.shape[0] // 2, :]
# 解析:式 (9.10) 换 |gamma| 即可(矩形线型的 gamma = sinc(.) e^{-j2pi nu0 tau},
# 等效于把余弦相位平移 2pi nu0 tau;对剖面形状只差整体平移,解析对比用同一约定)
I_an_r = analytic_two_beam(x2, 0.0, 650e-9, 0.25, 1e-3, 5e-3, deld, g_abs)
rel = np.max(np.abs(prof_r - I_an_r)) / I_an_r.max()
g_num = abs(np.sum(S_arr * np.exp(-2j * np.pi * nu_arr * tau)) * dnu)
W_sum = np.sum(S_arr) * dnu
print(f"谱线权重和 sum(S*dnu) = {W_sum:.4f}"
      f"(矩形边界落在样本之间,13 条线被矩形规则超计 -> > 1)")
print(f"离散度和 |sum S e^-j2pi nu tau dnu| = {g_num:.4f}(连续理论 {g_abs:.4f})")
print(f"仿真的条纹可见度 = 离散度和/权重和 = {g_num/W_sum:.4f}"
      f"(连续理论 {g_abs:.4f},差 {(g_num/W_sum-g_abs)/g_abs*100:+.1f}%)"
      f" —— 直接拿度和 0.8376 对比 0.827 会误读成 '+1.3%',先除权重和")
print(f"剖面 vs 式 (9.10)(|gamma| 换 {g_abs:.4f}):最大相对偏差 = {rel:.2e}")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].imshow(I_r, extent=[x2[0], x2[-1], x2[0], x2[-1]],
             origin="lower", cmap="gray")
ax[0].set_title("Ex 9.1: rect lineshape irradiance")
ax[1].plot(x2 * 1e4, prof_r)
ax[1].plot(x2[::3] * 1e4, I_an_r[::3], "o", ms=3, label=f"analytic |g|={g_abs:.3f}")
ax[1].set_xlim(-2, 2); ax[1].set_ylim(0, 1500)
ax[1].set_xlabel("x (1e-4 m)"); ax[1].legend(fontsize=8)
ax[1].set_title("Ex 9.1: profile vs analytic")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig9_ex91.png"), dpi=150)
OUT.append("out_fig9_ex91.png")

# ============================================================
# 习题 9.2:相干条纹位移(fringe shift)
# ============================================================
print("=" * 62)
print("习题 9.2 相干条纹位移")
M, L1 = 250, 50e-3
dx1 = L1 / M
x1, _ = coords(L1, M)
x1s = np.fft.fftshift(x1)
X1, Y1 = np.meshgrid(x1s, x1s)
w, dels, f = 1e-3, 5e-3, 0.25
lf = lam0 * f
_, fxc = coords(L1, M)
x2 = fxc * lf
fig, ax = plt.subplots(1, 4, figsize=(16, 3.2), sharey=True)
phases = []
per_ = f * lam0 / dels                            # 条纹周期 lam0 f / dels
# 包络基:像元化单束的 |G|^2(与仿真同一圆盘,自洽)
G1 = np.fft.fft2(circle(X1, Y1, w)) * dx1 ** 2 / lf
env = np.abs(np.fft.fftshift(G1))[M // 2, :] ** 2
win = np.abs(x2) < 2.5 * per_
for j, dd in enumerate([0, lam0 / 4, lam0 / 2, 3 * lam0 / 4]):
    u = circle(X1 - dels / 2, Y1, w) \
        + circle(X1 + dels / 2, Y1, w) * np.exp(1j * 2 * np.pi / lam0 * dd)
    u2 = np.fft.fft2(u) * dx1 ** 2 / lf
    I = np.fft.ifftshift(np.abs(u2) ** 2)
    pr = I[M // 2, :]
    # 条纹相位:prof = env*(a + b*cos(2pi x/per) + c*sin(2pi x/per)),
    # 相干情形该分解是精确的(包络与仿真同一像元化圆盘)
    th = 2 * np.pi * x2 / per_
    A3 = np.stack([env[win], env[win] * np.cos(th[win]),
                   env[win] * np.sin(th[win])], axis=1)
    co, *_ = np.linalg.lstsq(A3, pr[win], rcond=None)
    phases.append(np.arctan2(co[2], co[1]))
    ax[j].plot(x2 * 1e4, pr)
    ax[j].set_xlim(-1.5, 1.5)
    ax[j].set_title(f"$\\Delta d$ = {dd/lam0:.2f} $\\lambda_0$")
    ax[j].set_xlabel("x (1e-4 m)")
phases = np.array(phases)
frac = (-phases) / (2 * np.pi)          # 剖面 = cos(θ - φ_meas),理论相位 = +k0*dd
wrap = (frac + 0.5) % 1 - 0.5           # wrap 到 [-0.5, 0.5)
print(f"拟合条纹相位(周期数,wrap 到 [-0.5,0.5)):{np.round(wrap, 3)}")
print(f"理论 k0*dd/(2pi) mod 1:{np.round(np.array([0, .25, .5, -.25]), 3)}"
      f"(= deld/lam0 的 wrap;半周期翻转处平移取 ±0.5 等价)")
print("局限:图样以 lam0 为周期,deld 与 deld+lam0(整数倍波长)不可区分,"
      "单次测量只有 lam0 的无模糊范围;方向(超前/滞后)在一个周期内可分辨。")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig9_ex92.png"), dpi=150)
OUT.append("out_fig9_ex92.png")

# ============================================================
# 习题 9.3:可见度(式 9.31)+ 不等幅解析式 vs 仿真
# ============================================================
print("=" * 62)
print("习题 9.3 可见度 A1=0.75, A2=0.20")
A1, A2 = 0.75, 0.20
g_t = abs(gamma_temporal(tau, delnu, c / lam0))
V = 2 * A1 * A2 / (A1 ** 2 + A2 ** 2) * g_t
print(f"V = 2*A1*A2/(A1^2+A2^2) * |gamma(deld/c)| = "
      f"{2*A1*A2/(A1**2+A2**2):.4f} * {g_t:.4f} = {V:.4f}")
I_a, x2a, _ = pc_temp(amps=(A1, A2))
prof_a = I_a[I_a.shape[0] // 2, :]
# (c) 提议的解析式:式 (9.10) 的 [1+|g|cos] 括号换成
#     [A1^2+A2^2+2*A1*A2*|g|cos(...)];包络去掉 (9.10) 里双束等幅的因子 2
#     —— A1=A2=1 时括号 = 2[1+|g|cos],精确退回式 (9.10)。
#     幅度:|A|^2 = w^4/(lam0 f)^2 [J1(z)/(w~ r)]^2(单束夫琅禾费包络)
wtil = 1e-3 / (lam0 * 0.25)
r = np.abs(x2a)
zq = 2 * np.pi * wtil * r
from scipy.special import j1
ratio = np.where(r == 0, np.pi, j1(zq) / (wtil * np.where(r == 0, 1.0, r)))
env_a = 1 / (lam0 * 0.25) ** 2 * (1e-3 ** 2 * ratio) ** 2
I_an_a = env_a * (A1 ** 2 + A2 ** 2
                  + 2 * A1 * A2 * g_t * np.cos(2 * np.pi * 5e-3 * x2a / (lam0 * 0.25)
                                               + 2 * np.pi / lam0 * deld))
rel = np.max(np.abs(prof_a - I_an_a)) / I_an_a.max()
print(f"剖面 vs 提议解析式:最大相对偏差 = {rel:.2e}")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].imshow(I_a, extent=[x2a[0], x2a[-1], x2a[0], x2a[-1]],
             origin="lower", cmap="gray")
ax[0].set_title("Ex 9.3: A1=0.75, A2=0.20")
ax[1].plot(x2a * 1e4, prof_a)
ax[1].plot(x2a[::3] * 1e4, I_an_a[::3], "o", ms=3, label="proposed analytic")
ax[1].set_xlim(-2, 2); ax[1].set_ylim(0, 1500)
ax[1].set_xlabel("x (1e-4 m)"); ax[1].legend(fontsize=8)
ax[1].set_title(f"Ex 9.3: V = {V:.3f}")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig9_ex93.png"), dpi=150)
OUT.append("out_fig9_ex93.png")

# ============================================================
# 习题 9.4:l_cr = 15/8/4/1 mm 的 |mu_12| 与失效模式
# ============================================================
print("=" * 62)
print("习题 9.4 空间相干长度扫描")
lam, L1, M = 650e-9, 50e-3, 250
for Lcr in (15e-3, 8e-3, 4e-3, 1e-3):
    mu = np.exp(-(5e-3 / Lcr) ** 2)
    print(f"l_cr = {Lcr*1e3:4.0f} mm: |mu_12| = exp(-(5/l_cr)^2) = {mu:.4f}"
          f"(式 9.20 精确值, sigma_f=2.5 l_cr: {mu12_exact920(5e-3, Lcr, 2.5*Lcr):.4f})")
Nsc = 100
res = []
for Lcr, mult, n_use, tag in [
        (15e-3, 2.5, Nsc, "15mm, mult 2.5\n(sigma_f=37.5>L/2)"),
        (15e-3, 1.5, Nsc, "15mm, mult 1.5\n(sigma_f=22.5<L/2)"),
        (8e-3, 2.5, Nsc, "8mm, mult 2.5\n(book params)"),
        (4e-3, 2.5, Nsc, "4mm"),
        (1e-3, 2.5, Nsc, "1mm, N=100"),
        (1e-3, 2.5, 400, "1mm, N=400")]:
    Ij, xj, _ = pc_spatial(N=n_use, Lcr=Lcr, mult=mult, seed=7)
    Vj = fringe_visibility(Ij[Ij.shape[0] // 2, :], xj, 0.25 * lam / 5e-3)
    res.append((tag, Ij, xj, Vj))
    vs = "no fringes" if not np.isfinite(Vj) else f"{Vj:.3f}"
    print(f"  {tag.replace(chr(10), ' '):38s}: 可见度 = {vs}")
print("15mm/mult2.5:屏几乎恒定 -> 输出完全相干结果(书描述的失效模式);"
      "mult 1.5 后 sigma_f < L/2,可见度明显下降。"
      "1mm:l_cr << 束间距 5mm,|mu_12| = exp(-25) ~ 0 -> 无条纹, "
      "图样退化为双倍单束 Airy 斑;更多屏只让平均更平滑,不恢复条纹。")
fig, ax = plt.subplots(2, 6, figsize=(19, 6))
for j, (tag, Ij, xj, Vj) in enumerate(res):
    ext = [xj[0], xj[-1], xj[0], xj[-1]]
    ax[0, j].imshow(np.sqrt(Ij), origin="lower", extent=ext, cmap="gray")
    vs = "no fringes" if not np.isfinite(Vj) else f"V={Vj:.2f}"
    ax[0, j].set_title(tag + f"\n{vs}", fontsize=9)
    ax[1, j].plot(xj * 1e4, Ij[Ij.shape[0] // 2, :])
    ax[1, j].set_xlim(-2, 2)
    ax[1, j].set_xlabel("x (1e-4 m)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig9_ex94.png"), dpi=150)
OUT.append("out_fig9_ex94.png")

# ============================================================
# 习题 9.5:单个 5mm x 2mm 矩形孔,l_cr = 20 -> 1 mm
# ============================================================
print("=" * 62)
print("习题 9.5 矩形孔部分相干衍射")
fig, ax = plt.subplots(2, 4, figsize=(15, 6))
for j, Lcr in enumerate([20e-3, 8e-3, 4e-3, 1e-3]):
    Ij, xj, _ = pc_spatial(N=200, Lcr=Lcr, aperture="rect5x2", seed=7)
    ext = [xj[0], xj[-1], xj[0], xj[-1]]
    ax[0, j].imshow(np.sqrt(Ij), origin="lower", extent=ext, cmap="gray")
    ax[0, j].set_title(f"$\\ell_{{cr}}$ = {Lcr*1e3:.0f} mm")
    ax[1, j].plot(xj * 1e4, Ij[Ij.shape[0] // 2, :])
    ax[1, j].set_xlim(-2, 2)
    ax[1, j].set_xlabel("x (1e-4 m)")
print("解读:l_cr >> 孔径(20mm vs 5mm)时接近完全相干 -> sinc x sinc 图样带旁瓣;"
      "l_cr 减小,孔径内远距离点对的条纹贡献加权 exp(-(x'/l_cr)^2) 先死 -> "
      "旁瓣先消失、中央瓣变宽变滑;l_cr=1mm << 孔径 -> 图样趋于横向相干函数的"
      "傅里叶变换(光滑宽包络),孔径的衍射旁瓣完全消失。")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig9_ex95.png"), dpi=150)
OUT.append("out_fig9_ex95.png")

# ============================================================
# 习题 9.6:高斯谢尔模型光束(GSM)
# ============================================================
print("=" * 62)
print("习题 9.6 GSM 光束(w0=2cm, lam=1.06um, z=1km, l_cr=1cm, M=256)")
Mg, w0, lamg, zg, Lcrg = 256, 0.02, 1.06e-6, 1000.0, 0.01
L1g = 0.55                       # TF 判据 dx >= lam*z/L:2.148mm >= 1.927mm(书 5.4 节)
dxg = L1g / Mg
k = 2 * np.pi / lamg
xg, _ = coords(L1g, Mg)
# 注意:propTF(书 5.1 节)约定"源场按居中排列存放"(它内部自己做 fftshift),
# 所以这里用居中坐标生成 U1;相位屏 fie 是 ifft2 的输出(FFT 排列),
# 乘场之前要先 ifftshift 回居中——两套排列别混(见笔记坑与陷阱)。
Xg, Yg = np.meshgrid(xg, xg)
U1 = np.exp(-(Xg ** 2 + Yg ** 2) / w0 ** 2)          # 式 (9.34),居中排列
Delt = np.sqrt(1 + (2 * zg / (k * w0 ** 2)) ** 2
               * (1 + 2 * w0 ** 2 / Lcrg ** 2))      # 式 (9.33)
Delt_coh = np.sqrt(1 + (2 * zg / (k * w0 ** 2)) ** 2)
print(f"Delta = {Delt:.4f}(式 9.33;完全相干 {Delt_coh:.4f});"
      f"GSM 峰值 I0/Delta^2 = {1/Delt**2:.4f},相干 {1/Delt_coh**2:.4f}")
sf_g = 2.5 * Lcrg
sr_g = np.sqrt(4 * np.pi * sf_g ** 4 / Lcrg ** 2)
dfxg = 1 / L1g
_, fxg = coords(L1g, Mg)
FXg, FYg = np.meshgrid(np.fft.fftshift(fxg), np.fft.fftshift(fxg))
Fg = np.exp(-np.pi ** 2 * sf_g ** 2 * (FXg ** 2 + FYg ** 2))
Ng = 200
rng = np.random.default_rng(9)
I2g = np.zeros((Mg, Mg))
for n in range(Ng // 2):
    rt = rng.standard_normal((Mg, Mg)) + 1j * rng.standard_normal((Mg, Mg))
    fie = (np.fft.ifft2(Fg * rt) * sr_g / dfxg) * Mg ** 2 * dfxg ** 2
    fie_c = np.fft.ifftshift(fie)                    # FFT 排列 -> 居中排列
    u2 = propTF(U1 * np.exp(1j * fie_c.real), L1g, lamg, zg)
    I2g += np.abs(u2) ** 2
    u2 = propTF(U1 * np.exp(1j * fie_c.imag), L1g, lamg, zg)
    I2g += np.abs(u2) ** 2
I2g /= Ng
u2c = propTF(U1, L1g, lamg, zg)                      # 完全相干对比
I_coh = np.abs(u2c) ** 2
# 剖面与拟合
x2g = xg                                             # propTF 观察面与源面同格
pr_g = I2g[Mg // 2, :]
pr_c = I_coh[Mg // 2, :]
# 解析 GSM(式 9.32,I0=1:源面峰值辐照度,总功率 pi w0^2/2 守恒)
I_an_g = 1 / Delt ** 2 * np.exp(-2 * x2g ** 2 / (w0 * Delt) ** 2)
# 从仿真剖面拟合 Delta:ln I = ln(I0/D^2) - 2 r^2/(w0 D)^2
mfit = np.abs(x2g) < 0.12
sg = np.polyfit(x2g[mfit] ** 2, np.log(pr_g[mfit]), 1)[0]
D_fit = np.sqrt(2 / (-sg * w0 ** 2))
P_sim = np.sum(I2g) * dxg ** 2
print(f"拟合 Delta_sim = {D_fit:.3f}(理论 {Delt:.3f});"
      f"峰值:GSM 仿真 {pr_g.max():.4f} vs 解析 {1/Delt**2:.4f},"
      f"相干 {pr_c.max():.4f} vs {1/Delt_coh**2:.4f}")
print(f"总功率:GSM 仿真 {P_sim:.4e} m^2(源功率 pi*w0^2/2 = "
      f"{np.pi*w0**2/2:.4e},相位屏只调相不调幅,平均功率应守恒)")
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ext = [x2g[0], x2g[-1], x2g[0], x2g[-1]]
ax[0].imshow(np.sqrt(I2g), origin="lower", extent=ext, cmap="gray")
ax[0].set_title("Ex 9.6: GSM beam at z=1km (sqrt stretch)")
ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("y (m)")
ax[1].plot(x2g * 1e2, pr_g, label="GSM sim (N=200)")
ax[1].plot(x2g * 1e2, I_an_g, "--", label=f"GSM analytic, $\\Delta$={Delt:.3f}")
ax[1].plot(x2g * 1e2, pr_c, ":", label=f"coherent sim, $\\Delta$={Delt_coh:.3f}")
ax[1].plot(x2g * 1e2, np.exp(-2 * x2g ** 2 / (w0 * Delt_coh) ** 2) / Delt_coh ** 2,
           "-.", lw=0.8, alpha=0.6, label="coherent analytic")
ax[1].set_xlabel("x (cm)"); ax[1].legend(fontsize=8)
ax[1].set_title("Ex 9.6: profiles at z = 1 km")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "out_fig9_ex96.png"), dpi=150)
OUT.append("out_fig9_ex96.png")

print("=" * 62)
print("saved: " + ", ".join(OUT))
