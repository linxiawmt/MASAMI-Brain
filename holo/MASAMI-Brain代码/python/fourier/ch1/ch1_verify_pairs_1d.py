# ===== 书 1.5 节：用 numpy.fft 数值验证表 1.3 的一维变换对 =====
# 连续傅里叶变换的 DFT 近似（矩形求和）：
#     G(f) ≈ dx * Σ_n g(x_n) exp(-j2π f x_n)
# numpy.fft 的正变换不除 N（归一化在逆变换里），因此乘 dx 就回到连续 FT 的量纲。
# 信号定义在 x ∈ [-L/2, L/2) 上（0-based 索引），标准"居中"流程：
#     ifftshift -> fft -> fftshift
# 得到按频率从负到正排列的频谱（本书第 2、3 章系统讲这套定标）。
# MATLAB 注意：对偶数 N，fftshift 与 ifftshift 等价；MATLAB fft 同样不除 N。
import numpy as np
import matplotlib.pyplot as plt

N = 512            # 采样点数（偶数）
L = 8.0            # 空间窗口宽度 [m]
dx = L / N         # 采样间隔 [m]
n = np.arange(N) - N // 2
x = n * dx                       # x ∈ [-L/2, L/2)
f = n / (N * dx)                 # 频率轴 [cycles/m]，与 x 同序（居中）


def rect(t):
    """rect：|t|<1/2 为 1，|t|=1/2 为 1/2，其余 0（表 1.2）。"""
    t = np.asarray(t, dtype=float)
    r = np.where(np.abs(t) < 0.5, 1.0, 0.0)
    return np.where(np.abs(t) == 0.5, 0.5, r)


def tri(t):
    """Λ：|t|<=1 为 1-|t|，其余 0（表 1.2）。"""
    t = np.asarray(t, dtype=float)
    return np.where(np.abs(t) <= 1.0, 1.0 - np.abs(t), 0.0)


def ft_centered(g):
    """一维连续 FT 的 DFT 近似（居中输出）。"""
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(g))) * dx


def sincn(t):
    """归一化 sinc：sin(pi t)/(pi t)，与 numpy.sinc / MATLAB sinc 一致。"""
    return np.sinc(t)


# ---------- 对 1：rect(x) <-> sinc(f)（表 1.3 第 1 行，a=1） ----------
G1 = ft_centered(rect(x))
e1 = np.max(np.abs(G1 - sincn(f)))
print(f"[1] rect(x) <-> sinc(f)          最大绝对误差 = {e1:.3e}")

# ---------- 对 2：Λ(x) <-> sinc^2(f)（表 1.3 第 3 行，a=1） ----------
G2 = ft_centered(tri(x))
e2 = np.max(np.abs(G2 - sincn(f) ** 2))
print(f"[2] tri(x)  <-> sinc^2(f)        最大绝对误差 = {e2:.3e}")

# ---------- 对 3：exp(-pi x^2) <-> exp(-pi f^2)（表 1.3 第 5 行，a=1） ----------
g3 = np.exp(-np.pi * x**2)
G3 = ft_centered(g3)
e3 = np.max(np.abs(G3 - np.exp(-np.pi * f**2)))
print(f"[3] exp(-pi x^2) <-> exp(-pi f^2) 最大绝对误差 = {e3:.3e}")

# ---------- 定理 1：相似性 rect(x/2) <-> 2 sinc(2f)（表 1.1，a=2） ----------
G4 = ft_centered(rect(x / 2))
e4 = np.max(np.abs(G4 - 2 * sincn(2 * f)))
print(f"[4] 相似性 rect(x/2) <-> 2sinc(2f) 最大绝对误差 = {e4:.3e}")

# ---------- 定理 2：位移 rect(x-0.7) <-> sinc(f) exp(-j2*pi*f*0.7)（表 1.1） ----------
x0 = 0.7
G5 = ft_centered(rect(x - x0))
e5 = np.max(np.abs(G5 - sincn(f) * np.exp(-2j * np.pi * f * x0)))
print(f"[5] 位移   rect(x-0.7) 复数最大误差 = {e5:.3e}")

# ---------- 定理 3：中心纵坐标 G(0) = ∫g dx（表 1.1 最后一行） ----------
sums = np.sum(rect(x)) * dx
print(f"[6] 中心纵坐标  G(0): 数值 = {G1[N//2].real:.6f}, "
      f"求和*dx = {sums:.6f}, 解析 = 1.0")

# ---------- 勘误验证：表 1.3 最后一行 chirp 变换对（符号有误） ----------
# 原书：exp(-j*pi*x^2) <-> j*exp(-j*pi*f^2)（2D 版本 j*|ab|*exp[-j*pi(...)]）
# 正确：exp(-j*pi*x^2) <-> exp(-j*pi/4)*exp(+j*pi*f^2)
#       （2D 版本 -j*|ab|*exp[+j*pi(a^2 fX^2 + b^2 fY^2)]）
# chirp 不衰减，截断积分误差 ~ 1/(窗口半宽)，用宽窗口就能清晰区分两种公式。
Nc = 4096
Lc = 40.0
dxc = Lc / Nc
nc = np.arange(Nc) - Nc // 2
xc = nc * dxc
fc = nc / (Nc * dxc)
gc = np.exp(-1j * np.pi * xc**2)
Gc = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(gc))) * dxc
m = np.abs(fc) <= 1.0
mine = np.exp(-1j * np.pi / 4) * np.exp(1j * np.pi * fc**2)
book = 1j * np.exp(-1j * np.pi * fc**2)
print(f"[7] chirp 修正式误差 = {np.max(np.abs(Gc[m] - mine[m])):.3e}, "
      f"原书式误差 = {np.max(np.abs(Gc[m] - book[m])):.3e}  （应前者 << 后者）")

# ---------- 画对比图：数值（圆点） vs 解析（实线） ----------
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
for ax, G, A, ttl in [
    (axes[0], G1, sincn(f), r"rect $\leftrightarrow$ sinc"),
    (axes[1], G2, sincn(f) ** 2, r"tri $\leftrightarrow$ sinc$^2$"),
    (axes[2], G3, np.exp(-np.pi * f**2), r"gauss $\leftrightarrow$ gauss"),
]:
    m = np.abs(f) <= 4
    ax.plot(f[m], A[m], "b-", lw=1.5, label="analytic")
    ax.plot(f[m], G[m].real, "ro", ms=4, mfc="none", label="DFT numeric")
    ax.set_title(ttl); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    ax.set_xlabel("f [cycles/m]")
fig.suptitle(f"Table 1.3 pairs verified with numpy.fft (N={N}, L={L}, "
             f"max errors {e1:.1e} / {e2:.1e} / {e3:.1e})")
fig.tight_layout()
fig.savefig(__file__.replace("ch1_verify_pairs_1d.py", "out_1d_pairs.png"), dpi=150)
print("已保存 out_1d_pairs.png")
