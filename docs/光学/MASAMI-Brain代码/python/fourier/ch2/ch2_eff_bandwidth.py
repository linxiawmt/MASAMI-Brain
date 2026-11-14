# ===== 书 2.2 节：有效带宽（式 2.9-2.11 + 表 2.1） =====
# 有效带宽定义：包含 98% 总频谱功率的（径向）频谱宽度。
#   P(B)/P_T = 0.98,  P(B) = ∫_0^{2π}∫_0^B |F(ρcosθ, ρsinθ)|^2 ρ dρ dθ   （式 2.10）
# 总功率 P_T 用 Parseval 定理在空域算（式 2.9），数值累积一律除以解析 P_T，
# 避免频域截断把 98% 交点带偏（|F|^2 在 98% 附近非常平坦，0.3% 的归一化偏差
# 就能把 B 推移百分之几十）。
# 书中结果（表 2.1）：rect(x/2w)rect(y/2w) → B≈5/w；circ(sqrt(x^2+y^2)/w) → B≈5/w；
#                     exp[-π(x^2+y^2)/w^2] → B≈0.79/w。
# 圆函数的两种半径约定（承接 ch1 勘误 1）：
#   表 1.2 定义 circ(r/w) 是半径 w/2 的圆盘 → 径向 98% 半径 B·w ≈ 9.97；
#   Goodman/表 1.3 约定（半径 = w 的圆盘）→ B·w ≈ 4.99，与表 2.1 的 5/w 一致。
#   即表 2.1 的 circ 行再次站在"半径 = 分母"的 Goodman 约定一边。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import sinc, jinc   # 公共模块：rect/sinc/jinc 等，ch1 起的积木

w = 1.0                        # 特征半宽（结果按 B·w 报告，与 w 无关）
NR = 16001                     # ρ 网格密度
rho = np.linspace(0.0, 60.0 / w, NR)           # 频域半径 [1/m]


def cumtrapz(y, x):
    """0 起头的累积梯形积分（numpy 2.x 无内建 cumtrapz）。"""
    return np.concatenate(([0.0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(x))))


# ---------- rect：|F|^2 非圆对称，须做完整 θ 积分（式 2.10） ----------
# 利用正方形对称性：积分 θ∈[0, π/4] 再乘 8（被积函数在 π/4 镜像 + π/2 周期）。
NT = 16384
th = np.linspace(0.0, np.pi / 4, NT)
ring_rect = np.empty(NR)
step = 128
for s in range(0, NR, step):
    r = rho[s:s + step, None]
    fX = r * np.cos(th)[None, :]
    fY = r * np.sin(th)[None, :]
    v = (4.0 * w ** 2 * sinc(2 * w * fX) * sinc(2 * w * fY)) ** 2
    ring_rect[s:s + step] = (v[:, 1:] + v[:, :-1]).sum(axis=1) * (th[1] - th[0]) / 2
ring_rect *= 8.0
P_rect = cumtrapz(ring_rect * rho, rho)        # P(B)，除以解析 P_T=4w² 得占比

# ---------- circ 与 gauss：圆对称，θ 积分 trivially = 2π，直接径向 1D ----------
F2_circ = jinc(rho, w) ** 2                    # 半径 w 圆盘（Goodman/表 1.3 约定）
P_circ = cumtrapz(2 * np.pi * F2_circ * rho, rho)
F2_gauss = (w ** 2 * np.exp(-np.pi * w ** 2 * rho ** 2)) ** 2
P_gauss = cumtrapz(2 * np.pi * F2_gauss * rho, rho)

P_T = {"rect": 4.0 * w ** 2,                   # 式 2.9
       "circ(R=w)": np.pi * w ** 2,            # 圆盘面积
       "circ(R=w/2)": np.pi * (w / 2.0) ** 2,
       "gauss": w ** 2 / 2.0}                  # ∫∫exp(-2πr²/w²)

B_rect = np.interp(0.98 * P_T["rect"], P_rect, rho) * w
B_circ_w = np.interp(0.98 * P_T["circ(R=w)"], P_circ, rho) * w
B_circ_half = B_circ_w * 2                     # 半径减半 → B 加倍（尺度反比）
B_gauss = np.interp(0.98 * P_T["gauss"], P_gauss, rho) * w

# rect 交叉验证：直角网格 + 解析归一化（98% 附近的累积曲线极平坦，
# 极坐标的 θ/ρ 离散误差会放大成 B 的百分之几偏移；直角网格求和更稳）。
half, nc = 40.0, 4001
fc = np.linspace(-half, half, nc)
dfc = fc[1] - fc[0]
acc5 = acc55 = s_tot = 0.0
for s in range(0, nc, 250):
    FYc, FXc = np.meshgrid(fc[s:s + 250], fc, indexing="ij")
    Pc = 16.0 * sinc(2 * FXc) ** 2 * sinc(2 * FYc) ** 2
    Rc = np.hypot(FXc, FYc)
    s_tot += Pc.sum() * dfc ** 2
    acc5 += Pc[Rc <= 5.0].sum() * dfc ** 2
    acc55 += Pc[Rc <= 5.5].sum() * dfc ** 2
fr5, fr55 = acc5 / P_T["rect"], acc55 / P_T["rect"]
B_rect_cart = 5.0 + 0.5 * (0.98 - fr5) / (fr55 - fr5)   # 线性插值出 98% 交点

print("98% 有效带宽（数值，径向判据）")
print(f"  rect(x/2w)rect(y/2w) : B·w = {B_rect:.3f}（极坐标）；"
      f"直角网格参照 {B_rect_cart:.3f}  （书表 2.1：5/w ✓）")
print(f"  circ 半径 w   (Goodman/表 1.3 约定): B·w = {B_circ_w:.3f}（书表 2.1：5/w ✓）")
print(f"  circ 半径 w/2 (表 1.2 定义)        : B·w = {B_circ_half:.3f}")
print(f"  exp(-pi r^2/w^2)     : B·w = {B_gauss:.4f}  （书表 2.1：0.79/w；"
      f"解析 sqrt(ln50/2pi) = {np.sqrt(np.log(50) / (2 * np.pi)):.4f}）")

# circ 用 u=2πRρ 空间交叉验证：P/(πR²) = 2∫_0^{2πRρ} J1²(u)/u du（已知总积分 = 1/2），
# 一维精细梯形即可做到近机器精度（免二维截断误差）。
from scipy.special import j1 as _j1
u = np.linspace(0.0, 400.0, 400001)
vu = np.where(u > 1e-12, _j1(u) ** 2 / np.where(u > 1e-12, u, 1.0), 0.0)
cum_u = np.concatenate(([0.0], np.cumsum((vu[1:] + vu[:-1]) / 2 * np.diff(u))))
u98 = np.interp(0.98 * 0.5, cum_u, u)
print(f"  [交叉验证] u 空间精确解：98% 处 2πRρ = {u98:.2f} → B·R = {u98 / (2 * np.pi):.3f}"
      f"（半径 w 约定 B·w = {u98 / (2 * np.pi):.3f}；半径 w/2 约定 B·w = {u98 / np.pi:.3f}）")

# ---------- 画图：左两幅 = rect 频谱 98% 范围（复现图 2.3），右 = 累积功率曲线 ----------
fmax = 6.0 / w
n = 601
fx1 = np.linspace(-fmax, fmax, n)
FX, FY = np.meshgrid(fx1, fx1)                 # 'xy' 索引，与 MATLAB 一致
Fmag = 4.0 * w ** 2 * sinc(2 * w * FX) * sinc(2 * w * FY)

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
im = axes[0].imshow(np.abs(Fmag), extent=[-fmax, fmax, -fmax, fmax], origin="lower",
                    cmap="magma")
axes[0].add_patch(plt.Circle((0, 0), B_rect / w, fill=False, color="cyan", ls="--", lw=1.5))
axes[0].set_title(r"$|F(f_X,f_Y)|$, rect, circle = 98% power")
axes[0].set_xlabel(r"$f_X/w$"); axes[0].set_ylabel(r"$f_Y/w$")
plt.colorbar(im, ax=axes[0], fraction=0.046)

im = axes[1].imshow(np.abs(Fmag) ** 2, extent=[-fmax, fmax, -fmax, fmax], origin="lower",
                    cmap="magma")
axes[1].add_patch(plt.Circle((0, 0), B_rect / w, fill=False, color="cyan", ls="--", lw=1.5))
axes[1].set_title(r"$|F|^2$ (power spectrum), rect")
axes[1].set_xlabel(r"$f_X/w$"); axes[1].set_ylabel(r"$f_Y/w$")
plt.colorbar(im, ax=axes[1], fraction=0.046)

axes[2].plot(rho * w, P_rect / P_T["rect"], label="rect(x/2w)rect(y/2w)")
axes[2].plot(rho * w, P_circ / P_T["circ(R=w)"], label="circ, disk R=w (Table 2.1)")
axes[2].plot(rho * w, P_gauss / P_T["gauss"], label=r"$\exp(-\pi r^2/w^2)$")
axes[2].axhline(0.98, color="k", ls=":", lw=1)
axes[2].plot(B_rect, 0.98, "o", ms=5)
axes[2].plot(B_circ_w, 0.98, "o", ms=5)
axes[2].plot(B_gauss, 0.98, "o", ms=5)
axes[2].annotate(f"{B_rect:.2f}", (B_rect, 0.98), textcoords="offset points", xytext=(-6, 8), fontsize=8)
axes[2].annotate(f"{B_circ_w:.2f}", (B_circ_w, 0.98), textcoords="offset points", xytext=(8, -14), fontsize=8)
axes[2].annotate(f"{B_gauss:.2f}", (B_gauss, 0.98), textcoords="offset points", xytext=(4, -12), fontsize=8)
axes[2].set_xlim(0, 12); axes[2].set_ylim(0.9, 1.001)
axes[2].set_xlabel(r"$B\cdot w$"); axes[2].set_ylabel("cumulative power fraction")
axes[2].set_title("98% effective bandwidth (Table 2.1)")
axes[2].legend(fontsize=8); axes[2].grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_eff_bw.png"), dpi=150)
print("已保存 out_eff_bw.png")
