# -*- coding: utf-8 -*-
"""书 7.2.3 节 coh_image 相干成像例子的 Python 复现(书图 7.6/7.7/7.8)。

参数与书一致:wxp=6.25 mm、zxp=125 mm(f/10)、lambda=0.5 um、
L=0.3 mm、M=250(理想图像读 ch7_usaf_chart.py 生成的 USAF1951_py250.png,
替代原书网站已失效的 USAF1951B250.png)。
输出:out_fig7_6.png(相干传递函数)、out_fig7_7.png(模拟相干像)、
     out_fig7_8.png(v=-0.8e-4 m 行剖面)+ 命令行数值检查。

MATLAB/Python 差异:
  - MATLAB 偶数阵列中心样本在 (M/2+1,M/2+1)(1-based),numpy 在 (M//2,M//2);
    书第 49 行 vindex=round(vvalue/du+(M/2+1)) 对应本脚本 row。
  - 书第 30 行 surf(H.*.99) 的 *0.99 是 MATLAB 7.1 画单位高圆盘的显示
    权宜之计;matplotlib 无此问题。
  - 书第 34 行 xlabel('fv')/ylabel('fu') 与 surf(fu,fv,...) 的轴映射相反
    (surf 第一参数 fu 应配 x 轴)——圆瞳下无影响,矩形瞳(习题 7.3)
    必须改正,否则轴标反了。本脚本按正确映射标注。
  - imagesc + axis xy 对应 plt.imshow(origin='lower');
    nthroot(Ii,2) 对应 np.sqrt(Ii)。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
sys.path.insert(0, os.path.join(_here, "..", "common"))
from base import coords                                  # noqa: E402
from ch7_imaging import freq_coords, H_circ, coh_image   # noqa: E402
from ch7_usaf_chart import ensure_chart, L_IMG, N_PX     # noqa: E402

# ---------- 书第 3-12 行:读理想图像、归一化、坐标 ----------
A = ensure_chart() * 255.0            # 书 3 行 A=imread(...);0..255
A = np.flipud(A)                      # 书 5 行 A=flipud(A):底行 -> 负 v
Ig = A.astype(float)                  # 书 6-7 行:转浮点、峰值归一(已 0..255)
Ig = Ig / Ig.max()
ug = np.sqrt(Ig)                      # 书 9 行:辐照度图像 -> 场幅度(二值图无差别)

L = L_IMG                             # 书 10 行:0.3e-3 m
M = N_PX                              # 250
du = L / M                            # 书 11 行:1.2 um
u, _ = coords(L, M)                   # 书 12 行:u=-L/2:du:L/2-du
v = u

# ---------- 书第 19-27 行:成像系统参数与相干传递函数 ----------
lam = 0.5e-6                          # 书 19 行
wxp = 6.25e-3                         # 书 20 行:出瞳半径
zxp = 125e-3                          # 书 21 行:出瞳距离
f0 = wxp / (lam * zxp)                # 书 22 行:相干截止频率 1e5 cyc/m
fu = freq_coords(L, M)                # 书 24 行:频率坐标(间隔 1/L)
H = H_circ(fu, f0)                    # 书 27 行 H=circ(sqrt(Fu^2+Fv^2)/f0)

fN = 1.0 / (2 * du)                   # 理想图像的奈奎斯特频率

# ---------- 书第 35-39 行:式 (7.16) 相干成像 ----------
ui, Gi = coh_image(ug, H)
Ii = np.abs(ui) ** 2                  # 书 39 行

row = M // 2 + int(round(-0.8e-4 / du))   # 书 48-49 行:v=-0.8e-4 的行号

# ================= 出图(书图 7.6/7.7/7.8) =================
fig, ax = plt.subplots(figsize=(5.2, 4.4), subplot_kw={"projection": "3d"})
FU, FV = np.meshgrid(fu, fu)
ax.plot_surface(FU, FV, H * 0.99, cmap="gray", rstride=2, cstride=2,
                linewidth=0, antialiased=True)
ax.set_xlabel("fv (cyc/m)"); ax.set_ylabel("fu (cyc/m)")
ax.set_title("coherent transfer function H (f/10)")
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_6.png"), dpi=110)
plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.5))
axs[0].imshow(np.sqrt(Ii), extent=[u[0], u[-1] + du, v[0], v[-1] + du],
              origin="lower", cmap="gray", vmin=0, vmax=1)
axs[0].set_xlabel("u (m)"); axs[0].set_ylabel("v (m)")
axs[0].set_title("simulated coherent image (sqrt Ii)")
axs[0].annotate("Group -2\nElement 6", xy=(-0.8e-4, -0.8e-4),
                xytext=(-1.35e-4, -1.15e-4), color="w", fontsize=8,
                arrowprops=dict(arrowstyle="->", color="w"))
axs[1].plot(u, Ii[row, :], label="simulated")
axs[1].plot(u, Ig[row, :], ":", label="ideal")
axs[1].set_xlabel("u (m)"); axs[1].set_ylabel("Irradiance")
axs[1].set_title("row profile at v = -0.8e-4 m")
axs[1].legend(); axs[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_7.png"), dpi=110)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.6, 3.6))
ax.plot(u, Ii[row, :], label="simulated")
ax.plot(u, Ig[row, :], ":", label="ideal")
ax.set_xlabel("u (m)"); ax.set_ylabel("Irradiance")
ax.set_title("coh_image profile, v = -0.8e-4 m (Group -2 Element 6 / Element 1)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_8.png"), dpi=110)
plt.close(fig)

# ================= 数值检查(笔记 7.2.3 节) =================
print("== 书 7.2.3 节 coh_image 数值检查 ==")
print(f"du = {du*1e6:.2f} um,  判据 7.18: du <= lam*(f/#)/2 = {lam*10/2*1e6:.2f} um")
print(f"fN = 1/(2du)          = {fN:.3e} cyc/m")
print(f"f0 = wxp/(lam*zxp)    = {f0:.3e} cyc/m")
print(f"2f0 = 1/(lam*(f/#))   = {2*f0:.3e} cyc/m   (判据 7.17: 2f0 <= fN,裕量 fN/2f0 = {fN/(2*f0):.2f})")

# 场谱支撑:|Gi| > 1e-6*max 的最大半径(应 = f0,相干截止把谱拦住)
Gi_c = Gi                                   # 已转居中排列
fr = np.hypot(*np.meshgrid(fu, fu))
sup_f = fr[np.abs(Gi_c) > 1e-6 * np.abs(Gi_c).max()].max()
print(f"像面场谱支撑半径          = {sup_f:.3e} cyc/m  (理论 f0 = {f0:.3e})")

# 辐照度谱支撑:|F{Ii}| 的最大频率(应 = 2f0 —— 平方使频率内容翻倍)
FI = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Ii)))
sup_i = fr[np.abs(FI) > 1e-6 * np.abs(FI).max()].max()
print(f"像面辐照度谱支撑半径      = {sup_i:.3e} cyc/m  (理论 2f0 = {2*f0:.3e})")

# 组-2/单元6 条纹周期:剖面里 -1.1e-4 ~ -0.5e-4 m 区间的峰间距
seg = Ii[row, :]
uu = u
loc = (uu > -1.15e-4) & (uu < -0.45e-4)
pk = uu[loc][np.where(np.diff(np.sign(seg[loc] - 0.5 * seg[loc].max())) < 0)[0]]
if len(pk) >= 2:
    per = np.mean(np.diff(pk))
    print(f"组-2/单元6 条纹周期       = {per*1e6:.2f} um -> {1/per:.3e} cyc/m"
          f"  (书判读:周期约 1e-5 m,恰为相干截止 f0)")

# 振铃(Gibbs):组-2/单元1 大条(在 v=+1.0e-4 m 行)边缘附近的过冲
row_e1 = M // 2 + int(round(+1.0e-4 / du))
seg1 = Ii[row_e1, :][(uu > -1.14e-4) & (uu < -0.62e-4)]
print(f"组-2/单元1 边缘过冲峰值   = {seg1.max():.3f}(理想 1.0,相干振铃)")
print("已保存 out_fig7_6/7/8.png")
