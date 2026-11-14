# ===== 书 4.4.1/4.4.2 节：RS 传递函数 (4.19) vs 菲涅耳传递函数 (4.24) =====
# 对比两张"传递函数":看菲涅耳近似在哪个频率范围内可信、以及 |f|>1/lambda 的
# 倏逝波截止。相位差推导：
#   arg(H_RS/H_F) = kz(sqrt(1-(lam f)^2) - 1) + pi*lam*z*f^2
#     ≈ -pi*lam^3*z*f^4/4  （把式 (4.19) 的根号展开到 (lam f)^4 阶，二阶项与
#       式 (4.24) 的 -pi*lam*z*f^2 正好抵消——这是"二阶精度"近似的含义）
# MATLAB/Python 差异：np.sqrt 对负实数给出 nan 并 warning，要用复数输入
# (a + 0j) 才能"延拓"到倏逝区；MATLAB 的 sqrt(real) 同样只给实部/复数需 complex()。
import os

import numpy as np
import matplotlib.pyplot as plt

lam = 0.633e-6        # 波长 (m)，He-Ne
z = 50.0              # 传播距离 (m)（书 4.5 节例子的距离）
k = 2 * np.pi / lam

f = np.linspace(0, 1.5 / lam, 4000)      # 频率轴：以 1/lambda 为单位标定
prop = f < 1 / lam                        # 传播波分量（式 4.19 下方的条件）

# ---------- |H_RS|：传播区幅度 1，倏逝区指数衰减 ----------
H_RS = np.zeros_like(f, dtype=complex)
H_RS[prop] = np.exp(1j * k * z * np.sqrt(1 - (lam * f[prop]) ** 2 + 0j))
H_RS[~prop] = np.exp(-k * z * np.sqrt((lam * f[~prop]) ** 2 - 1))
mag_RS = np.abs(H_RS)

# ---------- 相位差（只定义在传播区） ----------
dphi = np.abs(k * z * (np.sqrt(1 - (lam * f[prop]) ** 2) - 1) + np.pi * lam * z * f[prop] ** 2)
f_asym = f[prop][f[prop] > 1e4]           # 渐近线 -pi*lam^3*z*f^4/4 只在大 f 画
asym = np.pi * lam ** 3 * z * f_asym ** 4 / 4
f_01 = (0.4 / (np.pi * lam ** 3 * z)) ** 0.25   # |dphi| = 0.1 rad 的频率

print(f"1/lambda（传播波截止）= {1 / lam:.3e} cyc/m")
print(f"|dphi|=0.1 rad 的频率 = {f_01:.3e} cyc/m = {f_01 * lam:.4f} /lambda")
print(f"例：f = 0.01/lambda 处相位差 = "
      f"{np.pi * lam ** 3 * z * (0.01 / lam) ** 4 / 4:.3f} rad"
      f"（f^4 增长：f 每大 10 倍，误差大 1 万倍）")

fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))

axes[0].semilogy(f * lam, mag_RS, label=r"$|H_{RS}|$  Eq. (4.19)")
axes[0].axhline(1.0, ls="--", lw=1, color="crimson", label=r"$|H_F| \equiv 1$  Eq. (4.24)")
axes[0].axvline(1.0, ls=":", lw=1, color="gray")
axes[0].set_ylim(1e-14, 3)
axes[0].set_xlabel(r"$f_X\ \lambda$"); axes[0].set_ylabel("magnitude")
axes[0].set_title(r"RS transfer function: evanescent cut at $|f|=1/\lambda$")
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)

axes[1].loglog(f[prop] * lam, dphi, label=r"$|\arg(H_{RS}/H_F)|$")
axes[1].loglog(f_asym * lam, asym, ":", lw=1, label=r"$\pi\lambda^3 z f_X^4/4$")
axes[1].axhline(0.1, ls="--", lw=1, color="crimson")
axes[1].axvline(f_01 * lam, ls=":", lw=1, color="gray")
axes[1].set_xlabel(r"$f_X\ \lambda$"); axes[1].set_ylabel("phase difference (rad)")
axes[1].set_title("Fresnel vs RS phase error (z = 50 m)")
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3, which="both")

fig.tight_layout()
fp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_transfer_compare.png")
fig.savefig(fp, dpi=150)
print("已保存 out_transfer_compare.png")
