# ===== 书 1.5 节：基本函数（复现原书图 1.1 / 图 1.2） =====
# 这些函数是全书建模"孔径 / 光场"的积木，定义必须与原书表 1.2 一致。
# MATLAB/Python 差异：numpy 索引从 0 开始（MATLAB 从 1 开始）；
# np.meshgrid 默认 indexing='xy' 与 MATLAB meshgrid 行为一致，
# 对圆对称函数（circ、高斯）两种索引结果相同。
import numpy as np
import matplotlib.pyplot as plt


# ---------- 一维基本函数（表 1.2） ----------
def rect(x):
    """rect(x)：|x| < 1/2 取 1，|x| = 1/2 取 1/2，其余取 0。"""
    x = np.asarray(x, dtype=float)
    r = np.where(np.abs(x) < 0.5, 1.0, 0.0)
    r = np.where(np.abs(x) == 0.5, 0.5, r)   # 边界点取半值；均匀网格很少恰好命中
    return r


def tri(x):
    """Λ(x)：|x| <= 1 取 1-|x|，其余取 0。"""
    x = np.asarray(x, dtype=float)
    return np.where(np.abs(x) <= 1.0, 1.0 - np.abs(x), 0.0)


def sinc(x):
    """sinc(x) = sin(pi*x)/(pi*x)；x=0 处极限为 1。
    注意 numpy 自带的 np.sinc 就是这个归一化定义，等价于本函数。"""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    mask = x != 0
    out[mask] = np.sin(np.pi * x[mask]) / (np.pi * x[mask])
    return out


def gauss(x):
    """高斯 exp(-pi*x^2)： Fourier 变换恰好是自己（表 1.3）。"""
    return np.exp(-np.pi * np.asarray(x, dtype=float) ** 2)


def circle(X, Y, R=0.5):
    """circ(sqrt(X^2+Y^2))：半径 R 的圆盘（表 1.2 中 R = 1/2）。
    注意：原书表 1.2 定义半径为 1/2，而表 1.3 的变换对对应半径 a 的圆盘，
    两者相差 2 倍半径——见 ch1_verify_circ_jinc.py 的数值验证。"""
    return np.where(np.hypot(X, Y) < R, 1.0, 0.0)


# ---------- 自检：定义是否与表 1.2 一致 ----------
print("rect(0)  =", rect(0.0), "(应为 1.0)")
print("rect(0.5)=", rect(0.5), "(应为 0.5)")
print("tri(0.5) =", tri(0.5), "(应为 0.5)")
print("sinc(0)  =", sinc(0.0), "(应为 1.0)")
print("sinc(0.5)=", sinc(0.5), "(应为 2/pi =", 2 / np.pi, ")")

# ---------- 复现原书图 1.1：一维基本函数 ----------
x = np.linspace(-1, 1, 801)
xw = np.linspace(-4, 4, 1201)

fig = plt.figure(figsize=(8, 9))
ax = fig.add_subplot(3, 2, 1)
ax.plot(x, rect(x)); ax.set_title("rect(x)"); ax.grid(True, alpha=0.3)
ax = fig.add_subplot(3, 2, 2)
ax.plot(xw, sinc(xw)); ax.set_title("sinc(x)"); ax.grid(True, alpha=0.3)
ax = fig.add_subplot(3, 2, 3)
ax.plot(x, tri(x)); ax.set_title(r"$\Lambda(x)$"); ax.grid(True, alpha=0.3)
ax = fig.add_subplot(3, 2, 4)
# comb 是 delta 周期阵列，网格上无法直接画，用 stem 箭头示意
xn = np.arange(-2, 3)
ax.stem(xn, np.ones_like(xn), basefmt="none")
ax.set_title("comb(x)"); ax.set_xlim(-2.5, 2.5); ax.set_ylim(0, 1.2)
ax = fig.add_subplot(3, 2, 5)
xg = np.linspace(-2, 2, 801)
ax.plot(xg, gauss(xg)); ax.set_title(r"$\exp(-\pi x^2)$"); ax.grid(True, alpha=0.3)
fig.suptitle("Basic 1D functions (reproduction of book Fig. 1.1)")
fig.tight_layout()
fig.savefig(__file__.replace("ch1_basic_functions.py", "out_1d_functions.png"), dpi=150)

# ---------- 复现原书图 1.2：二维 circ 与高斯 ----------
N = 80
xs = (np.arange(N) - N / 2) * (4.0 / N)
X, Y = np.meshgrid(xs, xs)               # 'xy' 索引，与 MATLAB meshgrid 一致
C = circle(X, Y, R=0.5)                  # 表 1.2：半径 1/2
G2 = np.exp(-np.pi * (X**2 + Y**2))      # 表 1.3 最后一类二维高斯（a=b=1）

fig = plt.figure(figsize=(10, 4.2))
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot_wireframe(X, Y, C, rstride=2, cstride=2, color="k", linewidth=0.4)
ax.set_title(r"$\mathrm{circ}(\sqrt{x^2+y^2})$")
ax = fig.add_subplot(1, 2, 2, projection="3d")
ax.plot_wireframe(X, Y, G2, rstride=2, cstride=2, color="k", linewidth=0.4)
ax.set_title(r"$\exp[-\pi(x^2+y^2)]$")
fig.suptitle("Examples of 2D functions (reproduction of book Fig. 1.2)")
fig.tight_layout()
fig.savefig(__file__.replace("ch1_basic_functions.py", "out_2d_functions.png"), dpi=150)
print("已保存 out_1d_functions.png / out_2d_functions.png")
