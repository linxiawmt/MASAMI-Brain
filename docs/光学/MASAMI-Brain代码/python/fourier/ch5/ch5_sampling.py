# -*- coding: utf-8 -*-
"""书 5.4 节 菲涅耳传播采样:artifacts 复现 + 判据数值验证 + 精度对拍。

内容(对应原书图 5.4/5.5/5.6 与表 5.1/5.2 的数值侧面):
  (1) 图 5.4 复现:TF/IR x 轴辐照度剖面,z = 1000-20000 m
      —— TF@4000 的 ticks、TF@20000 的 stair-step、IR@1000 的周期副本、
         IR@20000 的假旁瓣全部重现            -> out_fig5_4.png
  (2) 图 5.5 复现:对应的辐照度图样            -> out_fig5_5.png
  (3) 功率诊断(习题 5.3):源面/观察面功率比,采样 regime 的"免费体检表"
  (4) 表 5.1 判据的 B1 复核:方孔谱 98% 径向谱功率半径 vs 书值 5.98/w
  (5) 5.4.4 节精度对拍:式 (5.11) 方孔解析解(2500 点) vs 数值传播(250 点)
      —— log 幅度 + unwrap 相位,三案例:TF@1000 / TF@2000 / IR@20000
                                                   -> out_fig5_6.png
MATLAB/Python 差异:
  - scipy.special.fresnel(x) 返回 (S, C)——S 在前,与 MATLAB 符号工具箱的
    mfun('FresnelC'/'FresnelS') 各返回一个不同,搞反会把式 (5.11) 的实虚部对调。
  - MATLAB 的 unwrap 与 numpy 的 unwrap 行为一致(跳变 > pi 时加/减 2π 整数倍)。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import fresnel

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "common"))
from base import coords, sinc                       # 公共积木(ch1 起)
from ch5_prop import propTF, propIR

# ---------- 书 sqr_beam 参数(书 5.3 节) ----------
L = 0.5                                              # 边长 (m)
M = 250                                              # 样本数
dx = L / M                                           # 2 mm
x1, _ = coords(L, M)
lam = 0.5e-6                                         # 波长 (m)
k = 2 * np.pi / lam                                  # 波数
w = 0.051                                            # 源方孔半宽 (m)
X1, Y1 = np.meshgrid(x1, x1)
u1 = np.where((np.abs(X1) < w) & (np.abs(Y1) < w), 1.0, 0.0)   # 式 (5.9)
P1 = np.sum(np.abs(u1) ** 2) * dx ** 2               # 式 (5.25):源面功率

distances = [1000.0, 2000.0, 4000.0, 20000.0]        # 图 5.4/5.5 的四个距离
results = {}                                         # z -> (u2_tf, u2_ir)
for z in distances:
    results[z] = (propTF(u1, L, lam, z), propIR(u1, L, lam, z))

outdir = os.path.dirname(os.path.abspath(__file__))

# ---------- (1) 图 5.4 复现:剖面(左 TF 右 IR) ----------
fig, axs = plt.subplots(4, 2, figsize=(12.5, 13.0), sharex=True, sharey=True)
for i, z in enumerate(distances):
    for j, (tag, u2) in enumerate([("TF", results[z][0]), ("IR", results[z][1])]):
        ax = axs[i, j]
        ax.plot(x1, np.abs(u2[M // 2, :]) ** 2)
        ax.set_title(f"{tag} z = {z:g} m")
        ax.set_ylim(0, 2.1); ax.grid(alpha=0.3)
        if i == 3:
            ax.set_xlabel("x (m)")
        if j == 0:
            ax.set_ylabel("irradiance")
fig.tight_layout()
fig.savefig(os.path.join(outdir, "out_fig5_4.png"), dpi=110)
plt.close(fig)

# ---------- (2) 图 5.5 复现:辐照度图样 ----------
fig, axs = plt.subplots(4, 2, figsize=(9.0, 16.0))
for i, z in enumerate(distances):
    for j, (tag, u2) in enumerate([("TF", results[z][0]), ("IR", results[z][1])]):
        ax = axs[i, j]
        ax.imshow(np.abs(u2) ** 2,
                  extent=[x1[0], x1[-1] + dx, x1[0], x1[-1] + dx],
                  origin="lower", cmap="gray")
        ax.set_title(f"{tag} z = {z:g} m")
        ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(os.path.join(outdir, "out_fig5_5.png"), dpi=110)
plt.close(fig)

# ---------- (3) 功率诊断(习题 5.3):P2/P1 = sum(I2)dx^2 / P1 ----------
print("== 功率诊断(习题 5.3):P2/P1,源面功率 P1 = %.6f ==" % P1)
print("   z (m)    TF      IR")
for z in distances:
    p_tf = np.sum(np.abs(results[z][0]) ** 2) * dx ** 2 / P1
    p_ir = np.sum(np.abs(results[z][1]) ** 2) * dx ** 2 / P1
    print(f"  {z:7.0f}  {p_tf:.5f}  {p_ir:.5f}")

# ---------- (4) 表 5.1 判据的 B1 复核:方孔谱的径向谱功率占比 ----------
# 谱 |F|^2 = (4w^2)^2 sinc^2(2w fX) sinc^2(2w fY)(式 5.20);全功率(Parseval)
# = (2w)^2(式 2.9)。占比一律相对解析总功率归一化——sinc^2 尾部衰减慢,
# 若用"网格总和"当分母,截断误差会把 98% 交点带偏百分之几十(ch2 2.2 节同款教训)。
f1d = np.arange(-250.0, 250.0, 0.2)                  # cyc/m,零点间隔 1/(2w) ≈ 9.8
df = f1d[1] - f1d[0]
s1 = sinc(2 * w * f1d) ** 2                          # 一维 sinc^2(2w f)
G = np.outer(s1, s1) * (4 * w ** 2) ** 2 * df ** 2   # 二维可分离谱功率(dF^2 定标)
r = np.hypot(f1d[None, :], f1d[:, None])             # 半径网格
order = np.argsort(r.ravel())
cs = np.cumsum(G.ravel()[order]) / ((2 * w) ** 2)    # 按半径排序的累计占比(解析归一化)
r_s = r.ravel()[order]

def frac_in_disk(B):
    """半径 B 圆盘内的谱功率占比(相对解析总功率 (2w)^2)。"""
    i = np.searchsorted(r_s, B)
    return cs[i - 1]

B_book = 5.98 / w                                    # 式 (5.10):B1 = 5.98/w
B_ch2 = 5.0 / w                                      # ch2 表 2.1 / 式 (2.11):B ≈ 5/w
f98 = np.interp(0.98, cs, r_s)                       # 占比 = 98% 的半径
print("== B1 复核(方孔谱径向 98% 谱功率,解析归一化) ==")
print(f"占比 98% 的半径      B98      = {f98:.1f} cyc/m = {f98*w:.2f}/w (2wB = {2*w*f98:.2f})")
print(f"ch2 表 2.1 的 5/w = {B_ch2:.1f} cyc/m    占比 = {frac_in_disk(B_ch2)*100:.2f}%")
print(f"书式 (5.10) 5.98/w = {B_book:.1f} cyc/m  占比 = {frac_in_disk(B_book)*100:.2f}%")

# ---------- (5) 5.4.4 节精度对拍:式 (5.11) 解析 vs 数值 ----------
# 式 (5.11):U2 = e^{jkz}/(2j) * {DC + j DS}(x 部分) * {DC + j DS}(y 部分)
#   DC = C(xi1) - C(xi2),DS = S(xi1) - S(xi2)
#   xi1 = sqrt(2/(lam z))(-w - x), xi2 = sqrt(2/(lam z))(w - x)
# 注意 scipy 的 fresnel 返回 (S, C),S 在前!
cases = [(1000.0, "TF"), (2000.0, "TF"), (20000.0, "IR")]
fig, axs = plt.subplots(3, 2, figsize=(13.0, 14.0))
print("== 精度对拍(式 5.11,主体区 |x| <= L/4) ==")
print("   case            max 相位差(rad)   max 幅度相对差   unwrap 常数差(rad)")
for row_i, (z, method) in enumerate(cases):
    u2 = results[z][0] if method == "TF" else results[z][1]

    # 数值剖面(250 点)
    mag_n = np.abs(u2[M // 2, :])
    ph_n = np.unwrap(np.angle(u2[M // 2, :]))
    ph_n = ph_n - ph_n[M // 2]                       # 坑与陷阱:先减中心样本对齐

    # 解析剖面(2500 点)
    xa, dxa = coords(L, 2500)
    T = np.sqrt(2.0 / (lam * z))
    def _F(xc):                                      # 一维因子 [DC + j DS]
        xi1, xi2 = T * (-w - xc), T * (w - xc)
        S1, C1 = fresnel(xi1)
        S2, C2 = fresnel(xi2)
        return (C2 - C1) + 1j * (S2 - S1)
    U_ana = _F(xa) * _F(0.0) / (2j)                  # e^{jkz} 只差整体相位,丢弃
    mag_a = np.abs(U_ana)
    ph_a = np.unwrap(np.angle(U_ana))
    ph_a = ph_a - ph_a[len(xa) // 2]                 # 同样对齐中心

    const_shift = (np.angle(u2[M // 2, M // 2])
                   - np.angle(U_ana[len(xa) // 2]))  # 对齐前的中心常数差(信息量)
    mag_ni = np.interp(xa, x1, mag_n)                # 数值插值到解析细网格
    ph_ni = np.interp(xa, x1, ph_n)
    main = np.abs(xa) <= L / 4                       # 主体区
    d_mag = np.max(np.abs(mag_ni - mag_a)[main]) / np.max(mag_a)
    d_ph = np.max(np.abs(ph_ni - ph_a)[main])

    print(f"   {method} @ z={z:6.0f}    {d_ph:.4g}          {d_mag:.4g}        {const_shift:.2f}")

    # 图:左 log 幅度、右相位(复现原书图 5.6 布局)
    lab = f"{method} z={z:g} m (250 pts)"
    axs[row_i, 0].semilogy(xa, mag_a, "--", label="analytic (2500 pts)")
    axs[row_i, 0].semilogy(xa, mag_ni, "-", label=lab)
    axs[row_i, 0].set_title(f"({'ace'[row_i]}) {method} magnitude  z = {z:g} m")
    axs[row_i, 0].set_ylabel("Magnitude"); axs[row_i, 0].legend(fontsize=8)
    axs[row_i, 0].grid(alpha=0.3, which="both")
    axs[row_i, 1].plot(xa, ph_a, "--", label="analytic (2500 pts)")
    axs[row_i, 1].plot(xa, ph_ni, "-", label=lab)
    axs[row_i, 1].set_title(f"({'bdf'[row_i]}) {method} phase  z = {z:g} m")
    axs[row_i, 1].set_ylabel("Phase (rad)"); axs[row_i, 1].legend(fontsize=8)
    axs[row_i, 1].grid(alpha=0.3)
axs[2, 0].set_xlabel("x (m)"); axs[2, 1].set_xlabel("x (m)")
fig.tight_layout()
fig.savefig(os.path.join(outdir, "out_fig5_6.png"), dpi=110)
plt.close(fig)

print("已保存 out_fig5_4.png / out_fig5_5.png / out_fig5_6.png")
