# ===== 书 4.6 节习题 4.2-4.5 参考解答（数值部分） =====
# 4.2 三种衍射解的适用距离范围（解析估算 + 打印）；
# 4.3 验证菲涅耳传递函数 (4.24) = 菲涅耳脉冲响应 (4.23) 的傅里叶变换；
# 4.4 三个孔径的解析夫琅禾费场/辐照度：(a) 矩形孔 (b) 环形孔 (c) 双圆孔；
# 4.5 双圆孔两路振幅不同（A1=1, A2=0.4）的夫琅禾费辐照度。
# 解析式统一用书 4.4.3 节的写法：场 = e^{jkz}/(j lam z) * e^{jk(x^2+y^2)/(2z)}
#   * FT{U1}(fX=x/(lam z), fY=y/(lam z))；辐照度按式 (4.8) 取模方。
# MATLAB/Python 差异：MATLAB sinc 是工具箱函数（本目录脚本用 base.sinc 对应
# np.sinc 同定义）；besselj(1,u) 对应 scipy.special.j1(u)。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import coords, sinc, jinc   # 公共积木（ch1 起）


def jinc_book(u):
    """书 ch4 的 jinc：J1(2*pi*u)/u，u=0 -> pi（= base.jinc(u, R=1.0)）。"""
    return jinc(u, R=1.0)


lam = 0.633e-6   # 全部习题统一用 He-Ne 波长（习题 4.4 的指定值）

# ================= 习题 4.2：三种衍射解的适用距离 =================
# 孔径直径 1 mm -> 半径 w = 0.5 mm；习题 4.2 用 lam = 0.5 um。
w42 = 0.5e-3
lam42 = 0.5e-6
# 式 (4.27)：z^3 >> (pi/(4 lam)) * [(x-xi)^2+(y-eta)^2]^2_max
# 观察区取孔径的几何影区（直径 1 mm）-> x-xi、y-eta 各可达 1 mm
Dmax2 = (2 * w42) ** 2 + (2 * w42) ** 2
z_eq27 = (np.pi / (4 * lam42) * Dmax2 ** 2) ** (1 / 3)
z_nf1 = w42 ** 2 / lam42                    # NF = 1 的距离（Fresnel 区上界口径）
z_nf01 = 10 * w42 ** 2 / lam42              # NF = 0.1（书 4.5 节的工程判据）
z_fh = np.pi * (w42 ** 2 + w42 ** 2) / lam42   # 式 (4.29)：z >> k(xi^2+eta^2)_max/2
print("========== 习题 4.2（直径 1 mm 孔径，lam = 0.5 um） ==========")
print(f"RS：只需 r >> lam -> 任何宏观距离都满足（z > ~1 mm 即可）")
print(f"Fresnel：式(4.27)（观察区取几何影区）z >> {z_eq27 * 1e2:.1f} cm；"
      f"NF 视角：NF<1 对应 z > {z_nf1:.2f} m，平滑孔径可到 NF~20-30（z ~ {z_nf1 / 25:.3f} m）")
print(f"Fraunhofer：式(4.29) z >> {z_fh:.2f} m；工程口径 NF<0.1 -> z > {z_nf01:.1f} m")

# ================= 习题 4.3：H = FT{h} 的数值验证 =================
# h(x,y) = e^{jkz}/(j lam z) * exp[j k (x^2+y^2)/(2z)]   （式 4.23）
# H(fX,fY) = e^{jkz} * exp[-j pi lam z (fX^2+fY^2)]      （式 4.24）
# 采样要点：chirp 的瞬时频率 x/(lam z) 必须落在网格 Nyquist 内——这要求
# 比较频带 |f| <= L/(2 lam z)（驻点带），且它要低于 1/(2dx)（第 5 章判据预览）。
L3, M3, z3 = 0.02, 512, 2.0
dx3 = L3 / M3
x3, f3 = coords(L3, M3)                      # 空域/频域坐标（式 2.21/2.22）
X3, Y3 = np.meshgrid(x3, x3)
FX3, FY3 = np.meshgrid(f3, f3)
k3 = 2 * np.pi / lam
h = np.exp(1j * k3 * z3) / (1j * lam * z3) * np.exp(1j * k3 * (X3 ** 2 + Y3 ** 2) / (2 * z3))
H_num = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(h))) * dx3 ** 2   # FT ≈ dx^2 · DFT
H_an = np.exp(1j * k3 * z3) * np.exp(-1j * np.pi * lam * z3 * (FX3 ** 2 + FY3 ** 2))

band = L3 / (2 * lam * z3)                   # 驻点带 |f| <= L/(2 lam z)
inb8 = (np.abs(FX3) <= 0.8 * band) & (np.abs(FY3) <= 0.8 * band)
inb5 = (np.abs(FX3) <= 0.5 * band) & (np.abs(FY3) <= 0.5 * band)
dphi8 = np.angle(H_num[inb8] * np.conj(H_an[inb8]))
dphi5 = np.angle(H_num[inb5] * np.conj(H_an[inb5]))
print("\n========== 习题 4.3（L=0.02 m, M=512, z=2 m） ==========")
print(f"驻点带 |f| <= L/(2 lam z) = {band:.3e} cyc/m（Nyquist = {1 / (2 * dx3):.3e}）")
print(f"0.8 频带内最大缠绕相位差 = {np.abs(dphi8).max():.3f} rad")
print(f"0.5 频带内最大缠绕相位差 = {np.abs(dphi5).max():.3f} rad"
      f"（残差来自常数幅度 chirp 的窗口截断，第 5 章判据的预告）")

# ================= 习题 4.4：三个孔径的解析夫琅禾费图样 =================
# 各孔径参数与 z 的选取（式 4.29 的 ~10 倍裕量；L ≈ 5 倍中心特征宽度）：
#   (a) 矩形孔 wX=0.1 mm, wY=0.05 mm：I2 = (4wXwY/(lam z))^2 sinc^2 sinc^2
#   (b) 环形孔 w1=1 mm（外径）, w2=0.2 mm（遮拦）：
#         I2 = (w1^2/(lam z))^2 [jinc(w1 rho/(lam z)) - (w2/w1)^2 jinc(w2 rho/(lam z))]^2
#   (c) 双圆孔 w=1 mm, 中心距 Delta s=4 mm：
#         I2 = (2 w^2/(lam z))^2 cos^2(pi s x/(lam z)) jinc^2(w rho/(lam z))
cases = {
    "a": dict(wX=0.1e-3, wY=0.05e-3, z=5.0, L=0.16, M=512),
    "b": dict(w1=1e-3, w2=0.2e-3, z=50.0, L=0.2, M=512),
    "c": dict(w=1e-3, s=4e-3, z=500.0, L=0.6, M=512),
}

def grid(cs):
    x, _ = coords(cs["L"], cs["M"])
    X, Y = np.meshgrid(x, x)
    return x, X, Y

print("\n========== 习题 4.4 ==========")
# (a) 矩形孔
ca = cases["a"]
xa, Xa, Ya = grid(ca)
za, lza = ca["z"], lam * ca["z"]
Ia = (4 * ca["wX"] * ca["wY"] / lza) ** 2 * \
     sinc(2 * ca["wX"] * Xa / lza) ** 2 * sinc(2 * ca["wY"] * Ya / lza) ** 2
za_need = np.pi * (ca["wX"] ** 2 + ca["wY"] ** 2) / lam
row_a = ca["M"] // 2
zx = lza / (2 * ca["wX"])
print(f"(a) z=5 m：式(4.29)要求 z >> {za_need:.3f} m（裕量 {5 / za_need:.0f} 倍）；"
      f"峰 = (4wXwY/(lam z))^2 = {(4 * ca['wX'] * ca['wY'] / lza) ** 2:.3e}，"
      f"实测 {Ia[row_a, row_a]:.3e}")
print(f"    sinc 第一零点 x = lam z/(2wX) = {zx:.4f} m（y 向 {zx * ca['wX'] / ca['wY']:.4f} m）")

# (b) 环形孔
cb = cases["b"]
xb, Xb, Yb = grid(cb)
zb, lzb = cb["z"], lam * cb["z"]
rho_b = np.sqrt(Xb ** 2 + Yb ** 2)
amp_b = cb["w1"] ** 2 / lzb * (jinc_book(cb["w1"] * rho_b / lzb)
                              - (cb["w2"] / cb["w1"]) ** 2 * jinc_book(cb["w2"] * rho_b / lzb))
Ib = amp_b ** 2
zb_need = np.pi * (cb["w1"] ** 2 + cb["w2"] ** 2) / lam
row_b = cb["M"] // 2
prof_b = Ib[row_b, :]
peak_b_an = (np.pi * (cb["w1"] ** 2 - cb["w2"] ** 2) / lzb) ** 2
# 第一零点：解析 jinc 零点 0.61 lam z / w1 附近数值寻谷
i_b = int(np.argmin(np.abs(xb - 0.61 * lzb / cb["w1"])))
seg = prof_b[i_b - 8: i_b + 8]
print(f"(b) z=50 m：式(4.29)要求 z >> {zb_need:.2f} m（裕量 {50 / zb_need:.0f} 倍）；"
      f"NF = {cb['w1'] ** 2 / lzb:.4f}")
print(f"    峰 = (pi(w1^2-w2^2)/(lam z))^2 = {peak_b_an:.3e}，实测 {Ib[row_b, row_b]:.3e}"
      f"（无遮拦圆孔为 {(np.pi * cb['w1'] ** 2 / lzb) ** 2:.3e}，乘 (1-(w2/w1)^2)^2 = "
      f"{(1 - (cb['w2'] / cb['w1']) ** 2) ** 2:.3f}）")
print(f"    外环第一零点 ≈ 0.61*lam*z/w1 = {0.61 * lzb / cb['w1']:.5f} m，"
      f"谷值实测 {seg.min():.2e}")

# (c) 双圆孔
cc = cases["c"]
xc, Xc, Yc = grid(cc)
zc, lzc = cc["z"], lam * cc["z"]
rho_c = np.sqrt(Xc ** 2 + Yc ** 2)
env_c = cc["w"] ** 2 / lzc * jinc_book(cc["w"] * rho_c / lzc)   # 单孔包络（含 1/(lam z) 前置在下面）
Ic = (2 * env_c) ** 2 * np.cos(np.pi * cc["s"] * Xc / lzc) ** 2
zc_need = np.pi * ((cc["s"] / 2 + cc["w"]) ** 2 + cc["w"] ** 2) / lam
row_c = cc["M"] // 2
prof_c = Ic[row_c, :]
peak_c_an = (2 * np.pi * cc["w"] ** 2 / lzc) ** 2
print(f"(c) z=500 m：式(4.29)要求 z >> {zc_need:.1f} m（裕量 {500 / zc_need:.0f} 倍）；"
      f"NF = {cc['w'] ** 2 / lzc:.4f}")
print(f"    峰 = (2 pi w^2/(lam z))^2 = {peak_c_an:.3e}，实测 {Ic[row_c, row_c]:.3e}")
print(f"    包络（单孔 Airy）第一零点 0.61*lam*z/w = {0.61 * lzc / cc['w']:.4f} m；"
      f"条纹周期 lam*z/Delta s = {lzc / cc['s']:.5f} m")

# ================= 习题 4.5：两路振幅不同（A1=1, A2=0.4） =================
# U2 = (e^{jkz}/(j lam z)) e^{chirp} * (w^2/(lam z)) jinc(w rho/(lam z))
#      * (A1 e^{-j pi s x/(lam z)} + A2 e^{+j pi s x/(lam z)})
# 提示（书原话）：先取模方，再把两个复指数合并成余弦：
#   I2 = (w^2/(lam z))^2 jinc^2(...) * (A1^2 + A2^2 + 2 A1 A2 cos(2 pi s x/(lam z)))
A1, A2 = 1.0, 0.4
Ie = (cc["w"] ** 2 / lzc) ** 2 * jinc_book(cc["w"] * rho_c / lzc) ** 2 * \
     (A1 ** 2 + A2 ** 2 + 2 * A1 * A2 * np.cos(2 * np.pi * cc["s"] * Xc / lzc))
prof_e = Ie[row_c, :]
vis = 2 * A1 * A2 / (A1 ** 2 + A2 ** 2)
print(f"\n========== 习题 4.5（A1=1, A2=0.4, z=500 m） ==========")
print(f"峰值 = (w^2/(lam z))^2 pi^2 (A1+A2)^2 = "
      f"{(cc['w'] ** 2 / lzc) ** 2 * np.pi ** 2 * (A1 + A2) ** 2:.3e}，实测 {Ic[row_c, row_c] * 0 + prof_e.max():.3e}")
print(f"条纹可见度 V = 2A1A2/(A1^2+A2^2) = {vis:.4f}；"
      f"中心区 Imax/Imin = ((A1+A2)/(A1-A2))^2 = {((A1 + A2) / (A1 - A2)) ** 2:.2f}")
print(f"条纹周期 lam*z/Delta s = {lzc / cc['s']:.5f} m")

# ================= 画图 =================
fig1, axes1 = plt.subplots(1, 2, figsize=(10.5, 4.2))
# 4.3：缠绕相位差剖面 + 低频区两相位叠加
row3 = M3 // 2
dphi_line = np.angle(H_num[row3, :] * np.conj(H_an[row3, :]))
sel = np.abs(f3) <= 0.8 * band
axes1[0].plot(f3[sel], dphi_line[sel])
axes1[0].set_xlabel(r"$f_X$ (cyc/m)"); axes1[0].set_ylabel("wrapped phase diff (rad)")
axes1[0].set_title(r"Ex 4.3: $\angle(H_{num}/H_{an})$ along $f_X$ ($f_Y=0$)")
axes1[0].grid(True, alpha=0.3)
sel2 = np.abs(f3) <= 400
axes1[1].plot(f3[sel2], np.unwrap(np.angle(H_num[row3, sel2])), label=r"$\angle$ FFT$\{h\}$")
axes1[1].plot(f3[sel2], np.unwrap(np.angle(H_an[row3, sel2])), ":", lw=2,
              label=r"$\angle H$  Eq.(4.24)")
axes1[1].set_xlabel(r"$f_X$ (cyc/m)"); axes1[1].set_ylabel("unwrapped phase (rad)")
axes1[1].set_title("low-freq window: the two chirps coincide")
axes1[1].legend(fontsize=8)
axes1[1].grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_ex43.png"), dpi=150)

fig2, axes2 = plt.subplots(2, 3, figsize=(12.5, 7.2))
for ax, I, cs, ttl in [
        (axes2[0, 0], Ia, ca, "(a) rect: $I_2^{1/3}$ image"),
        (axes2[1, 0], Ib, cb, "(b) annulus: $I_2^{1/3}$ image"),
        (axes2[1, 1], Ic, cc, "(c) pair: $I_2^{1/3}$ image")]:
    xg, _ = coords(cs["L"], cs["M"])
    ax.imshow(I ** (1 / 3), extent=[xg[0], xg[-1], xg[0], xg[-1]],
              origin="lower", cmap="gray")
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_title(ttl)

axes2[0, 1].plot(xa, Ia[row_a, :])
axes2[0, 1].set_xlabel("x (m)"); axes2[0, 1].set_ylabel("Irradiance")
axes2[0, 1].set_title("(a) profile: $(4w_Xw_Y/\\lambda z)^2\\,\\mathrm{sinc}^2\\mathrm{sinc}^2$")
axes2[0, 1].grid(True, alpha=0.3)

axes2[0, 2].plot(xb, prof_b)
axes2[0, 2].set_xlabel("x (m)"); axes2[0, 2].set_ylabel("Irradiance")
axes2[0, 2].set_title("(b) profile: annulus vs unobscured")
env1 = (cb["w1"] ** 2 / lzb) ** 2 * jinc_book(cb["w1"] * xb / lzb) ** 2
axes2[0, 2].plot(xb, env1, ":", lw=1)
axes2[0, 2].legend(["annulus", "unobscured disk"], fontsize=8)
axes2[0, 2].grid(True, alpha=0.3)

axes2[1, 2].plot(xc, prof_c)
axes2[1, 2].set_xlabel("x (m)"); axes2[1, 2].set_ylabel("Irradiance")
axes2[1, 2].set_title(r"(c) profile: fringes $\times$ Airy envelope")
axes2[1, 2].grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_ex44.png"), dpi=150)

fig3, axes3 = plt.subplots(1, 2, figsize=(10.5, 4.2))
xg, _ = coords(cc["L"], cc["M"])
axes3[0].imshow(Ie ** (1 / 3), extent=[xg[0], xg[-1], xg[0], xg[-1]],
                origin="lower", cmap="gray")
axes3[0].set_xlabel("x (m)"); axes3[0].set_ylabel("y (m)")
axes3[0].set_title(r"Ex 4.5: $A_1=1, A_2=0.4$ — $I_2^{1/3}$ image")
axes3[1].plot(xc, prof_e)
axes3[1].set_xlabel("x (m)"); axes3[1].set_ylabel("Irradiance")
axes3[1].set_title(r"profile: visibility $2A_1A_2/(A_1^2+A_2^2)=0.69$")
axes3[1].grid(True, alpha=0.3)
fig3.tight_layout()
fig3.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_ex45.png"), dpi=150)

print("\n已保存 out_ex43.png, out_ex44.png, out_ex45.png")
