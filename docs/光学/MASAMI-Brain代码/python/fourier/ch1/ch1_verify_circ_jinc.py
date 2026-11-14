# ===== 书 1.5 节：circ 的二维傅里叶变换对验证（原书图 1.3 的 "jinc" 曲面） =====
# 表 1.2 定义 circ(sqrt(x^2+y^2)) 在半径 1/2 内取 1；
# 表 1.3 却给出 circ(r/a) <-> a^2*J1(2*pi*a*rho)/(a*rho)，峰值 pi*a^2（a=1 时为 pi，
# 正是图 1.3 标注的峰值）。半径 a/2 的圆盘面积是 pi*a^2/4，变换峰值应等于面积
# （中心纵坐标定理），所以表 1.3 这一对应的是"半径 a 的圆盘"（Goodman 约定）。
# 下面用二维 DFT 同时验证两种半径，看解析式到底匹配哪一个。
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1

N = 512            # 每边采样点数
L = 4.0            # 空间窗口宽度 [m]
dx = L / N
n = np.arange(N) - N // 2
x = n * dx
f = n / (N * dx)                 # 频率轴 [cycles/m]
X, Y = np.meshgrid(x, x)         # 'xy' 索引（与 MATLAB meshgrid 一致）；圆对称下无所谓


def disk(R):
    """半径 R 的圆盘（边界测度为零，取值不影响数值积分）。"""
    return np.where(np.hypot(X, Y) < R, 1.0, 0.0)


def ft2_centered(g):
    """二维连续 FT 的 DFT 近似（居中输出）：G ≈ dx^2 * Σ g exp(-j2π(f·r))。"""
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(g))) * dx**2


def jinc(rho, R):
    """半径 R 圆盘的解析频谱：R*J1(2*pi*R*rho)/rho，rho=0 处极限 = 圆盘面积。"""
    out = np.empty_like(rho)
    nz = rho != 0
    out[nz] = R * j1(2 * np.pi * R * rho[nz]) / rho[nz]
    out[~nz] = np.pi * R**2
    return out


results = {}
for R in (1.0, 0.5):
    G = ft2_centered(disk(R))
    row = G[N // 2, :].real       # fY = 0 的径向剖面（圆对称，行/列等价）
    A = jinc(f, R)
    err = np.max(np.abs(row - A))
    area = np.sum(disk(R)) * dx**2          # 中心纵坐标定理：G(0,0) = 面积
    results[R] = (f, row, A, err)
    print(f"R = {R:>4}: 频谱峰值(数值) = {row[N//2]:.6f}, 解析极限 = {np.pi*R**2:.6f}, "
          f"圆盘面积(数值) = {area:.6f}, 径向剖面最大误差 = {err:.3e}")

print("\n结论：解析式 R*J1(2*pi*R*rho)/rho 精确匹配'半径 R 的圆盘'。")
print("表 1.3（a=1，峰值 pi）对应半径 1 的圆盘；若坚持表 1.2 的半径 1/2 定义，")
print("变换应为 0.5*J1(pi*rho)/rho（峰值 pi/4）——两处约定不一致，使用时务必小心。")

# ---------- 画径向剖面对比图 ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, R, ttl in [
    (axes[0], 1.0, "R = 1  (peak $\\pi$, matches Table 1.3 / Fig. 1.3)"),
    (axes[1], 0.5, "R = 1/2 (peak $\\pi/4$, matches Table 1.2 definition)"),
]:
    f_, row, A, err = results[R]
    m = np.abs(f_) <= 3.0
    ax.plot(f_[m], A[m], "b-", lw=1.5, label="analytic  $R\\,J_1(2\\pi R\\rho)/\\rho$")
    ax.plot(f_[m], row[m], "ro", ms=4, mfc="none", label="2D DFT numeric")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_title(f"{ttl}\nmax error = {err:.1e}")
    ax.set_xlabel(r"$\rho$ [cycles/m]"); ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
fig.suptitle("Circle function transform (jinc) — radial cut fY = 0")
fig.tight_layout()
fig.savefig(__file__.replace("ch1_verify_circ_jinc.py", "out_circ_jinc.png"), dpi=150)
print("已保存 out_circ_jinc.png")
