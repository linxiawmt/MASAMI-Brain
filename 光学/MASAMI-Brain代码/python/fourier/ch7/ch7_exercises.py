# -*- coding: utf-8 -*-
"""书 7.4 节习题 7.1-7.10 的 Python 实跑(解析题直接算,仿真题出图)。

7.1/7.2/7.8      解析:透镜定律、放大率、f/#、采样判据(命令行输出)。
7.3              矩形出瞳(wx=1 mm、wy=0.5 mm):x/y 方向截止不同 ->
                 相干/非相干像的各向异性分辨(谱支撑 + 图像)。
7.4              环形出瞳(wo=1 mm、wi=0.5 mm):相干传递函数挖掉 DC 与
                 低频 -> 边缘增强;OTF 与全孔径对比。
7.6              衍射限 PSF:h = F^-1{circ(rho/f0)} = pi f0^2 2J1(2 pi f0 u)/(2 pi f0 u),
                 PSF = |h|^2 为 Airy 图样;第一零点 = 1.22 lam (f/#)。
7.7              Rayleigh 判据:delta = 1.22 lam (f/#) = 12.2 um,
                 双点源间隔 S = 10 个样本(12 um),看凹陷。
7.9              相衬成像:beta = pi/2 正相衬(Ii ~ 1+2*theta)、
                 3*pi/2 负相衬(1-2*theta)、其他 beta。
7.10             傅里叶全息:记录 |FFT{物体+参考点}|^2,回放 IFFT,
                 挖掉轴上亮点;像 + 共轭像 + 中间自相关。

输出:out_fig7_ex73/74/7677/79/710.png + 命令行数值检查。

MATLAB/Python 差异:
  - 习题 7.9 的 DC 点:MATLAB G(1,1)(1-based),numpy G[0,0](0-based)。
  - 习题 7.10 的 udelta 来自附录 C.6(本目录 ch7_imaging.py 提供);
    MATLAB 侧 udelta.m 同源。
  - numpy 翻转 180 度写作 a[::-1, ::-1],MATLAB 是 rot90(a,2)。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import j1

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
sys.path.insert(0, os.path.join(_here, "..", "common"))
from base import coords, rect, ft2_centered              # noqa: E402  公共积木
from ch7_imaging import (freq_coords, H_circ, H_rect, H_annulus,   # noqa: E402
                         otf_from_H, coh_image, incoh_image,
                         otf_circ_analytic, udelta)
from ch7_usaf_chart import ensure_chart, L_IMG, N_PX     # noqa: E402

lam = 0.5e-6

print("== 习题 7.1(解析)==")
f = 50e-3; z1 = 200e-3
z2 = 1.0 / (1.0 / f - 1.0 / z1)
Mt = -z2 / z1
DXP = 2 * 10e-3
print(f"(a) z2 = {z2*1e3:.2f} mm, Mt = {Mt:.4f}(倒像)")
print(f"(b) 薄透镜:DXP = {DXP*1e3:.0f} mm, zXP = z2 = {z2*1e3:.2f} mm, "
      f"f/# = {z2/DXP:.3f}")
print(f"(c) z1 >> f: 1/z2 = 1/f - 1/z1 -> z2 -> {f*1e3:.0f} mm(物体趋远,像聚到焦面)")

print("== 习题 7.2(解析)==")
DXP2 = 0.5 * 0.0254; zXP2 = 4 * 0.0254; lam2 = 0.587e-6
fn2 = zXP2 / DXP2
f0_2 = 1.0 / (2 * lam2 * fn2)
du2 = lam2 * fn2 / 2
print(f"(a) f/# = {fn2:.1f}, 相干截止 f0 = wXP/(lam zXP) = {f0_2:.4e} cyc/m, "
      f"非相干 2f0 = {2*f0_2:.4e} cyc/m")
print(f"(b) du <= lam*(f/#)/2 = {du2*1e6:.3f} um")
print(f"(c) L = 1 mm -> M >= L/du = {1e-3/du2:.1f} -> 至少 426 样本(实用取 512)")

print("== 习题 7.8(解析,物空间 f/#)==")
fOBJ = 700e3 / 1.0
two_f0 = 1.0 / (lam * fOBJ)
delta8 = 1.22 * lam * fOBJ
du8 = lam * fOBJ / 2
print(f"(a) f/#OBJ = zEP/DEP = {fOBJ:.2e}")
print(f"(b) 非相干截止 2f0 = 1/(lam f/#OBJ) = {two_f0:.4f} cyc/m; "
      f"Rayleigh delta = 1.22*lam*(f/#OBJ) = {delta8:.3f} m")
print(f"(c) du <= lam*(f/#OBJ)/2 = {du8:.3f} m; 1024^2 -> L <= {1024*du8:.1f} m")

# ================= 习题 7.3:矩形出瞳 =================
L3, M3 = 1.0e-3, 250
du3 = L3 / M3
u3, _ = coords(L3, M3)
fu3 = freq_coords(L3, M3)
wx3, wy3, zxp3 = 1.0e-3, 0.5e-3, 50e-3
f0x = wx3 / (lam * zxp3)          # 4e4 cyc/m
f0y = wy3 / (lam * zxp3)          # 2e4 cyc/m
fN3 = 1.0 / (2 * du3)

Ig = ensure_chart(); Ig = np.flipud(Ig); Ig /= Ig.max()

H3 = H_rect(fu3, fu3, f0x, f0y)            # 相干传递函数(方通带)
ui3, _ = coh_image(np.sqrt(Ig), H3)
Ii3c = np.abs(ui3) ** 2
OTF3 = otf_from_H(H3)                      # 矩形瞳 OTF = Lambda*Lambda(式 7.29 推广)
Ii3i = incoh_image(Ig, OTF3)

FI3 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Ii3i)))
col_max = np.abs(FI3).max(axis=0)          # 每个 fu 列的最大幅度
row_max = np.abs(FI3).max(axis=1)
sup_x = np.abs(fu3)[col_max > 1e-6 * col_max.max()].max()
sup_y = np.abs(fu3)[row_max > 1e-6 * row_max.max()].max()

fig = plt.figure(figsize=(11.0, 8.6))
ax1 = fig.add_subplot(2, 2, 1, projection="3d")
FU3, FV3 = np.meshgrid(fu3, fu3)
ax1.plot_surface(FU3, FV3, H3 * 0.99, cmap="gray", rstride=2, cstride=2,
                 linewidth=0, antialiased=True)
ax1.set_title(f"coherent TF: f0x={f0x:.1e}, f0y={f0y:.1e} cyc/m")
ax2 = fig.add_subplot(2, 2, 2, projection="3d")
ax2.plot_surface(FU3, FV3, np.fft.fftshift(np.abs(OTF3)), cmap="gray",
                 rstride=2, cstride=2, linewidth=0, antialiased=True)
ax2.set_title("OTF (rect pupil)")
ax3 = fig.add_subplot(2, 2, 3)
ax3.imshow(np.sqrt(Ii3c), extent=[u3[0], u3[-1] + du3, u3[0], u3[-1] + du3],
           origin="lower", cmap="gray")
ax3.set_title("coherent image (rect pupil)")
ax4 = fig.add_subplot(2, 2, 4)
ax4.imshow(np.sqrt(Ii3i), extent=[u3[0], u3[-1] + du3, u3[0], u3[-1] + du3],
           origin="lower", cmap="gray")
ax4.set_title("incoherent image (rect pupil)")
for ax in (ax3, ax4):
    ax.set_xlabel("u (m)"); ax.set_ylabel("v (m)")
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_ex73.png"), dpi=105)
plt.close(fig)

print("== 习题 7.3 数值检查 ==")
print(f"du = {du3*1e6:.1f} um; 判据 7.18: u 方向 du <= lam*(zXP/2wx)/2 = "
      f"{lam*zxp3/(2*wx3)/2*1e6:.1f} um,v 方向 {lam*zxp3/(2*wy3)/2*1e6:.1f} um -> 均满足")
print(f"(b)(c) 相干像:|Fu| <= {f0x:.1e}、|Fv| <= {f0y:.1e} 的方通带")
print(f"非相干像谱支撑:u 方向 {sup_x:.3e}(理论 2f0x = {2*f0x:.3e}), "
      f"v 方向 {sup_y:.3e}(理论 2f0y = {2*f0y:.3e})")
print(f"(d) 本板条纹(映射到 L=1 mm):组-1/单元2 = 3.8e4 cyc/m < f0x -> u 向可分辨;"
      f"组-2/单元2 = 1.9e4 < f0y=2e4 但组-2/单元3 = 2.1e4 > f0y -> v 向只到组-2/单元2"
      f" ——竖条清楚、横向细节被抹掉,各向异性确凿")

# ================= 习题 7.4:环形出瞳 =================
wo3, wi3 = 1.0e-3, 0.5e-3
f0o = wo3 / (lam * zxp3)                   # 4e4(外半径决定相干截止)
f0i = wi3 / (lam * zxp3)                   # 2e4(内半径决定挖掉的 DC 半径)
H4 = H_annulus(fu3, f0o, f0i)
ui4, _ = coh_image(np.sqrt(Ig), H4)
Ii4c = np.abs(ui4) ** 2
ui4f, _ = coh_image(np.sqrt(Ig), H_circ(fu3, f0o))   # 同外径全孔径对照
OTF4 = otf_from_H(H4)
Ii4i = incoh_image(Ig, OTF4)

bar = np.zeros((M3, M3), bool)
bar[194:232, 32:72] = True                 # 组-2/单元1 三条竖条区域(图板行 18-55 翻转)
amp_full = np.abs(ui4f)[bar].mean()
amp_ann = np.abs(ui4)[bar].mean()

otf_ann_c = np.fft.fftshift(np.real(OTF4))[M3 // 2, :]
otf_full_c = otf_circ_analytic(fu3, f0o)[M3 // 2, :]
sl4 = np.abs(fu3) < 1.2e5
mid = (np.abs(fu3) > 2.5e4) & (np.abs(fu3) < 3.5e4)
high = (np.abs(fu3) > 5.5e4) & (np.abs(fu3) < 6.5e4)

fig = plt.figure(figsize=(11.0, 8.6))
ax1 = fig.add_subplot(2, 2, 1, projection="3d")
ax1.plot_surface(FU3, FV3, H4 * 0.99, cmap="gray", rstride=2, cstride=2,
                 linewidth=0, antialiased=True)
ax1.set_title(f"coherent TF annulus: f0o={f0o:.1e}, f0i={f0i:.1e} cyc/m")
ax2 = fig.add_subplot(2, 2, 2, projection="3d")
ax2.plot_surface(FU3, FV3, np.fft.fftshift(np.abs(OTF4)), cmap="gray",
                 rstride=2, cstride=2, linewidth=0, antialiased=True)
ax2.set_title("OTF (annular pupil)")
ax3 = fig.add_subplot(2, 2, 3)
ax3.imshow(np.sqrt(Ii4c), extent=[u3[0], u3[-1] + du3, u3[0], u3[-1] + du3],
           origin="lower", cmap="gray")
ax3.set_title("coherent image (annular: DC+low freq removed)")
ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(fu3[sl4], otf_ann_c[sl4], label="annular OTF (fV=0 slice)")
ax4.plot(fu3[sl4], otf_full_c[sl4], "--", label="full aperture OTF")
ax4.set_xlabel("fu (cyc/m)"); ax4.legend(); ax4.grid(alpha=0.3)
ax4.set_title("OTF slice: annular vs full")
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_ex74.png"), dpi=105)
plt.close(fig)

print("== 习题 7.4 数值检查 ==")
print(f"相干截止由外半径决定:f0 = {f0o:.3e}(内半径挖掉 |f| < {f0i:.3e} 的 DC 与低频)")
print(f"(d) 大条内部平均场:全孔径 {amp_full:.3f} -> 环孔 {amp_ann:.3f}"
      f"(平坦区域的场被抽掉,只剩边缘贡献 -> 边缘增强/伪相衬外观)")
print(f"OTF 环孔/全孔:低频(|f|~3e4)= "
      f"{otf_ann_c[mid].mean()/otf_full_c[mid].mean():.2f},"
      f"中频(|f|~6e4)= {otf_ann_c[high].mean()/otf_full_c[high].mean():.2f}"
      f"(低频受损、中高频相对占优)")

# ================= 习题 7.6:衍射限 PSF(Airy)=================
L6, M6 = L_IMG, N_PX
du6 = L6 / M6
u6, _ = coords(L6, M6)
f06 = 1.0e5                                          # f/10 系统
H6 = H_circ(freq_coords(L6, M6), f06)
df6 = 1.0 / L6                                       # H 的频率采样间隔
h6 = ft2_centered(H6, df6)                           # H 实偶 -> FT = IFT(base 公共积木)
row6 = M6 // 2
prof_h = np.real(h6[row6, :])
h0_theory = np.pi * f06 ** 2
h_analytic = np.pi * f06 ** 2 * 2 * j1(2 * np.pi * f06 * u6) / (2 * np.pi * f06 * u6)
h_analytic[M6 // 2] = h0_theory

PSF6 = np.abs(h6) ** 2
prof_psf = PSF6[row6, :] / PSF6[row6, M6 // 2]
airy = (2 * j1(2 * np.pi * f06 * u6) / (2 * np.pi * f06 * u6)) ** 2
airy[M6 // 2] = 1.0
hline = np.real(h6[row6, M6 // 2:])                  # h 沿 +u
nz = int(np.argmax(hline < 0))                       # 第一过零(样本)
# 线性插值细化过零位置:零点在 nz-1 与 nz 号样本之间
first_zero = ((nz - 1) + hline[nz - 1] / (hline[nz - 1] - hline[nz])) * du6

print("== 习题 7.6 数值检查 ==")
print(f"h(0) 数值 = {np.real(h6[row6, M6//2]):.4e}, 理论 pi*f0^2 = {h0_theory:.4e}")
print(f"(b)(c) PSF 第一零点 = {first_zero*1e6:.2f} um, 理论 1.22*lam*(f/#) = "
      f"{1.22*lam*10*1e6:.2f} um -> 中央瓣全宽 2.44*lam*(f/#) = {2*first_zero*1e6:.2f} um")
print(f"PSF vs Airy [2J1(v)/v]^2 最大差 = {np.max(np.abs(prof_psf - airy)):.2e}(归一化后)")

# ================= 习题 7.7:Rayleigh 判据 =================
fn7 = 20.0
zxp7 = 125e-3
wxp7 = zxp7 / (2 * fn7)                              # 3.125 mm
f07 = wxp7 / (lam * zxp7)                            # 5e4
delta7 = 1.22 * lam * fn7                            # 12.2 um
S7 = int(round(delta7 / du6))                        # 10.17 -> 10
Ig7 = np.zeros((M6, M6))
Ig7[M6 // 2, M6 // 2 - S7 // 2] = 1.0
Ig7[M6 // 2, M6 // 2 + S7 // 2] = 1.0
OTF7 = otf_from_H(H_circ(freq_coords(L6, M6), f07))
Ii7 = incoh_image(Ig7, OTF7)
prof7 = Ii7[M6 // 2, :]
pk = max(prof7[M6 // 2 - S7 // 2], prof7[M6 // 2 + S7 // 2])
dip_ratio = prof7[M6 // 2] / pk

fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.6))
axs[0].plot(u6, prof_h / h0_theory, label="numeric h(u)")
axs[0].plot(u6, h_analytic / h0_theory, "--", label="analytic 2J1(v)/v")
axs[0].set_title("Ex 7.6: impulse response h(u)")
axs[0].set_xlabel("u (m)"); axs[0].legend(); axs[0].grid(alpha=0.3)
axs[0].set_xlim(-3e-5, 3e-5)
axs[1].plot(u6, prof_psf, label="numeric PSF")
axs[1].plot(u6, airy, "--", label="Airy [2J1(v)/v]^2")
axs[1].set_title("Ex 7.6: PSF (Airy pattern)")
axs[1].set_xlabel("u (m)"); axs[1].legend(); axs[1].grid(alpha=0.3)
axs[1].set_xlim(-3e-5, 3e-5)
axs[2].plot(u6, prof7)
axs[2].set_title(f"Ex 7.7: two points, S={S7} px = {S7*du6*1e6:.0f} um, "
                 f"dip = {dip_ratio:.2f}")
axs[2].set_xlabel("u (m)"); axs[2].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_ex7677.png"), dpi=105)
plt.close(fig)

print("== 习题 7.7 数值检查 ==")
print(f"(a) f/# = {fn7:g} -> wxp = {wxp7*1e3:.3f} mm, delta = 1.22*lam*(f/#) = "
      f"{delta7*1e6:.2f} um")
print(f"(b) S = delta/du = {delta7/du6:.2f} -> 取 {S7} 个样本(间隔 {S7*du6*1e6:.0f} um)")
print(f"(d) 中央凹陷/峰 = {dip_ratio:.3f}(Airy 理论 0.735 -> 恰好'刚可分辨')")
print("已保存 out_fig7_ex7677.png")

# ================= 习题 7.9:相衬成像 =================
theta = np.sqrt(Ig) * np.pi / 100.0                  # 书 (a):theta = sqrt(Ig)*pi/100
ug9 = np.exp(1j * theta)
Ii_direct = np.abs(ug9) ** 2                         # 书 (b):直接算辐照度


def phase_contrast(beta):
    """书 (c)/(d):FFT -> DC 点乘 exp(j*beta) -> IFFT(保持图像方向)。"""
    G = np.fft.fft2(ug9)                             # DC 在 [0,0](numpy 0-based)
    G[0, 0] = G[0, 0] * np.exp(1j * beta)
    return np.abs(np.fft.ifft2(G)) ** 2


Ii_pos = phase_contrast(np.pi / 2)
Ii_neg = phase_contrast(3 * np.pi / 2)
Ii_pi = phase_contrast(np.pi)

row9 = M6 // 2
fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.6))
axs[0, 0].imshow(Ii_direct, extent=[u6[0], u6[-1] + du6, u6[0], u6[-1] + du6],
                 origin="lower", cmap="gray", vmin=1 - 1e-6, vmax=1 + 1e-6)
axs[0, 0].set_title(f"(b) direct |ug|^2 (dev from 1: {np.abs(Ii_direct-1).max():.1e})")
axs[0, 1].imshow(np.sqrt(Ii_pos), extent=[u6[0], u6[-1] + du6, u6[0], u6[-1] + du6],
                 origin="lower", cmap="gray")
axs[0, 1].set_title(r"(c) positive contrast, $\beta=\pi/2$")
axs[1, 0].imshow(np.sqrt(Ii_neg), extent=[u6[0], u6[-1] + du6, u6[0], u6[-1] + du6],
                 origin="lower", cmap="gray")
axs[1, 0].set_title(r"(d) negative contrast, $\beta=3\pi/2$")
seg = (u6 > -1.14e-4) & (u6 < -0.62e-4)              # 组-2 大条区域剖面
axs[1, 1].plot(u6[seg], Ii_pos[row9, seg], label=r"$\beta=\pi/2$")
axs[1, 1].plot(u6[seg], Ii_neg[row9, seg], label=r"$\beta=3\pi/2$")
axs[1, 1].axhline(1.0, color="k", lw=0.5)
axs[1, 1].plot(u6[seg], 1 + 2 * theta[row9, seg], ":", label=r"ideal $1\pm2\theta$")
axs[1, 1].set_title("profiles across Group -2 bars")
axs[1, 1].set_xlabel("u (m)"); axs[1, 1].legend(); axs[1, 1].grid(alpha=0.3)
for ax in (axs[0, 0], axs[0, 1], axs[1, 0]):
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_ex79.png"), dpi=105)
plt.close(fig)

on_pos = Ii_pos[row9, seg][theta[row9, seg] > 0].mean()
off_pos = Ii_pos[row9, seg][theta[row9, seg] == 0].mean()
on_neg = Ii_neg[row9, seg][theta[row9, seg] > 0].mean()
off_neg = Ii_neg[row9, seg][theta[row9, seg] == 0].mean()

print("== 习题 7.9 数值检查 ==")
print(f"(b) 直接辐照度 = 1(最大偏离 {np.abs(Ii_direct-1).max():.1e},浮点级)")
print(f"(c) beta=pi/2 正相衬:条纹上 {on_pos:.4f} vs 条纹外 {off_pos:.4f},"
      f"差 = {on_pos-off_pos:+.4f}(理论 2*theta = +{2*np.pi/100:.4f})")
print(f"(d) beta=3pi/2 负相衬:条纹上 {on_neg:.4f} vs 条纹外 {off_neg:.4f},"
      f"差 = {on_neg-off_neg:+.4f}(理论 -2*theta)")
print(f"(e) beta=pi:DC 反相 -> 线性衬度消失:max|Ii-1| = {np.abs(Ii_pi-1).max():.2e}"
      f"(只剩二阶小量;beta=0 时 Ii 恒为 1)")
print("已保存 out_fig7_ex79.png")

# ================= 习题 7.10:傅里叶全息 =================
M0, L0 = 500, 10e-3
dx0 = L0 / M0                                        # 0.02 mm
x0v, _ = coords(L0, M0)                              # -5 mm : 0.02 mm : ...
X0, Y0 = np.meshgrid(x0v, x0v)
mm = 1e-3
object10 = (rect((X0 + 2.0 * mm) / (0.4 * mm)) * rect(Y0 / (2.0 * mm))    # 竖梁
            + rect((X0 + 1.3 * mm) / (1.0 * mm)) * rect((Y0 - 0.8 * mm) / (0.4 * mm))
            + rect((X0 + 1.3 * mm) / (1.0 * mm)) * rect(Y0 / (0.4 * mm)))
point10 = 1000.0 * udelta((X0 - 2.0 * mm) / dx0) * udelta(Y0 / dx0)
u1_10 = object10 + point10                           # 书 (a)(b)

G10 = np.fft.fft2(u1_10)                             # 书 (c):记录透镜 = FFT
T10 = np.abs(G10) ** 2                               # 胶片记录的辐照度
P10 = np.fft.ifft2(T10)                              # 书 (d):回放(逆 FFT 保方向)
P10[0, 0] = 0.0                                      # 书 (e):挖掉轴上亮点(DC)

xc = x0v * 1e3                                       # mm(显示用)


def ncc(a, b):
    """归一化互相关系数(窗口匹配度)。"""
    a = a - a.mean(); b = b - b.mean()
    return float(np.sum(a * b) / np.sqrt(np.sum(a * a) * np.sum(b * b)))


r_w = int(round(2.2e-3 / dx0 / 2))                   # 2.2 mm 半宽 -> 55 px
Pabs = np.fft.fftshift(np.abs(P10))                  # 回放场(居中排列)

# 物体自身窗口(中心 x=-1.5 mm、y=0)与 180 度翻转版
cx_obj = int(round(-1.5e-3 / dx0)) + M0 // 2
obj_c = object10[M0 // 2 - r_w:M0 // 2 + r_w, cx_obj - r_w:cx_obj + r_w]
flip_c = obj_c[::-1, ::-1]


def window10(center_x_mm):
    c = int(round(center_x_mm * 1e-3 / dx0)) + M0 // 2
    return Pabs[M0 // 2 - r_w:M0 // 2 + r_w, c - r_w:c + r_w]


wL = window10(-3.5)                                  # 真像(物中心 -1.5 - x_ref)
wR = window10(+3.5)                                  # 共轭像(翻转)

fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.6))
film_log = np.log10(np.abs(np.fft.fftshift(G10)) ** 2
                    / np.abs(G10).max() ** 2 + 1e-12)
axs[0].imshow(film_log, vmin=-8, vmax=0,
              extent=[xc[0], xc[-1], xc[0], xc[-1]], cmap="gray")
axs[0].set_title("(c) film: log|FFT{object+ref}|^2")
axs[1].imshow(np.log10(Pabs + 1e-1),
              extent=[xc[0], xc[-1], xc[0], xc[-1]], cmap="gray")
axs[1].set_title("(f) playback log|IFFT{T}| (DC removed)")
axs[2].imshow(np.sqrt(Pabs / Pabs.max()),
              extent=[xc[0], xc[-1], xc[0], xc[-1]], cmap="gray")
axs[2].set_title("playback, sqrt scale (image / conjugate / autocorr)")
for ax in axs:
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_ex710.png"), dpi=105)
plt.close(fig)

print("== 习题 7.10 数值检查 ==")
print(f"物体(F 形,中心 -1.5 mm)+ 参考点(1000 倍,(2,0) mm)已建,图像 {M0}x{M0}")
print(f"左侧窗口(-3.5 mm)与物体 NCC = {ncc(wL, obj_c):+.3f},与翻转物体 = "
      f"{ncc(wL, flip_c):+.3f} -> 真像")
print(f"右侧窗口(+3.5 mm)与物体 NCC = {ncc(wR, obj_c):+.3f},与翻转物体 = "
      f"{ncc(wR, flip_c):+.3f} -> 共轭像")
print("中间 = 物体自相关(u1 ⊛ u1,傅里叶自相关定理);轴上亮点 = R^2*delta 已清零")
print("已保存 out_fig7_ex73/74/710.png 全部完成")
