# ===== 书 3.2-3.5 节：fft_example 脚本（建向量 -> 移位 -> FFT 定标 -> 解析对比） =====
# 书把这几节写成一个逐步加长的脚本 fft_example.m（原书第 1-58 行），本文件
# 按同样的参数与步骤复现书图 3.2-3.8，并附数值检查。
# 参数（书 3.2 节）：矩形半宽 w = 0.055 m（全宽 0.11 m）、边长 L = 2 m、
# 样本数 M = 200 -> 采样间隔 dx = 0.01 m。
# MATLAB/Python 差异：
#   - numpy 索引 0-based：居中排列的零坐标在 index M//2（本例 100），
#     MATLAB 是 index M/2+1（本例 101）；
#   - 书用 fftshift 做 FFT 前移位（偶数 M 下与 ifftshift 等价）；本文件按
#     第 2 章的严谨配对写 ifftshift -> fft -> fftshift，奇数 M 也不出错；
#   - MATLAB 的 plot '-o'、axis([...]) 对应 plt.plot(..., 'o-') 与
#     plt.xlim/plt.ylim。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import rect, coords, sinc   # 公共积木：表 1.2 函数与网格约定（ch2 起）


def rect_book(x):
    """书 3.1 节的 rect 编码：abs(x) <= 1/2 -> 1（边界取 1，不取半值）。

    表 1.2 的解析定义（base.rect）在边界取半值 1/2；两者只在恰好命中
    边界的样本上不同。书 3.1 节原话：半值会被解读成"斜坡"，不如单点
    跳变；代价是矩形全宽总由奇数个样本组成（请求偶数宽会多给一个样本），
    所以全书例子故意用 0.055 这类非整半宽。
    """
    return (np.abs(x) <= 0.5).astype(float)


# ---------- 3.2 创建向量 ----------
w = 0.055          # 矩形半宽 (m)。故意不取 0.05：见"奇数样本"特性
L = 2.0            # 向量边长 (m)
M = 200            # 样本数
dx = L / M         # 采样间隔 (m)
x, fx = coords(L, M)           # 式 (2.21)/(2.22)：x ∈ [-1, 0.99]，fx ∈ [-50, 49.5]
f = rect_book(x / (2 * w))     # 采样矩形：|x| <= 0.055 的样本为 1（11 个）

print(f"x 范围: [{x[0]:.2f}, {x[-1]:.2f}]，共 {len(x)} 个样本")
print(f"f==1 的样本数 = {int(f.sum())}（全宽 0.11 m / dx 0.01 m = 11，奇数）")
print(f"坐标 x=0 的 index = {int(np.argmin(np.abs(x)))}"
      f"（Python 0-based；MATLAB 是 {M // 2 + 1} = M/2+1）")

# 两种 rect 约定对照：本例网格上完全相同（边界 x=±0.055 不在网格上）
diff = np.max(np.abs(f - rect(x / (2 * w))))
print(f"rect_book 与表 1.2 的 base.rect 在本网格最大差异 = {diff:.1f}"
      f"（边界样本不在网格 -> 相同；若请求偶数全宽就会多出一个样本）")

# ---------- 3.3 为 FFT 移位 ----------
f0 = np.fft.ifftshift(f)   # 居中 -> 移位排列（零坐标样本挪到 index 0）
# 移位后中心样本在 index 0，5 个样本在"右"（index 1..5），其余 5 个在
# "左"——挂在向量末尾（index 195..199），对应书图 3.5。
# 不移位直接 fft 的后果：相当于对平移了 -L/2 的函数做变换，频域多出一线性
# 相位因子 exp(+j*2*pi*f*L/2)（位移定理）——幅度不变、相位全错。

# ---------- 3.4 计算 FFT 并定标 ----------
F0 = np.fft.fft(f0) * dx   # FFT 并乘 dx 回连续量纲（式 2.19 + 黎曼和因子）
F = np.fft.fftshift(F0)    # 居中排列
# fx 已由 coords 给出（式 2.22）：-50 : 0.5 : 49.5

# ---------- 3.5 与解析结果对比 ----------
F_an = 2 * w * sinc(2 * w * fx)   # rect(x/2w) <-> 2w sinc(2w fx)（表 1.3，a = 2w）
# np.sinc 与书用的 MATLAB 内置 sinc（信号处理工具箱）同定义：sin(pi u)/(pi u)

# ---------- 数值检查 ----------
ic = int(np.argmax(np.abs(F)))
print(f"DC 样本 |F(fx=0)| = {abs(F[ic]):.6f}（= 矩形面积 2w = {2 * w:.2f}，中心纵坐标定理）")
in_main = np.abs(fx) <= 25
err_all = np.max(np.abs(np.abs(F) - np.abs(F_an)))
err_main = np.max(np.abs(np.abs(F[in_main]) - np.abs(F_an[in_main])))
print(f"|FFT|-|解析| 最大差异：全频带 {err_all:.4f}；|fx|<=25 内 {err_main:.2e}")
print("  （边缘差异 = 周期延拓带来的混叠，书 3.5 节原话'FFT 边缘略高'）")
print(f"sinc 零点间隔 1/(2w) 折合样本数：L/(2w) = {L / (2 * w):.2f}（书说'约每 18 个点'）")
print(f"峰值处相位 = {np.angle(F[ic]):.2e}（实偶函数 -> 相位只有 0 与 ±π，±π 同一相位 mod 2π）")

# ---------- 画图：复现书图 3.3-3.8 ----------
fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.8))

axes[0, 0].plot(x, f, "o-", ms=4, color="crimson", mfc="none")     # 书图 3.3
axes[0, 0].set_xlim(-0.2, 0.2); axes[0, 0].set_ylim(0, 1.5)
axes[0, 0].set_xlabel("x (m)")
axes[0, 0].set_title("f = rect(x/(2w)) samples (Fig. 3.3)")

axes[0, 1].plot(np.arange(1, M + 1), f, "o-", ms=4, mfc="none")    # 书图 3.4（1-based 展示）
axes[0, 1].set_xlim(80, 120); axes[0, 1].set_ylim(0, 1.5)
axes[0, 1].set_xlabel("index (MATLAB 1-based)")
axes[0, 1].set_title("vs index: center at 101 = M/2+1 (Fig. 3.4)")

axes[0, 2].plot(np.arange(1, M + 1), f0, "o-", ms=4, mfc="none")   # 书图 3.5
axes[0, 2].set_xlim(0, 200); axes[0, 2].set_ylim(0, 1.5)
axes[0, 2].set_xlabel("index (MATLAB 1-based)")
axes[0, 2].set_title("after shift: center at index 1 (Fig. 3.5)")

axes[1, 0].plot(np.abs(F0))                                        # 书图 3.6(a)
axes[1, 0].set_xlabel("index")
axes[1, 0].set_title(r"$|F_0|$ vs index (Fig. 3.6a)")

axes[1, 1].plot(np.angle(F0))                                      # 书图 3.6(b)
axes[1, 1].set_xlabel("index")
axes[1, 1].set_title("phase vs index: only 0 and $\\pm\\pi$ (Fig. 3.6b)")

axes[1, 2].plot(fx, np.abs(F), fx, np.abs(F_an), ":")              # 书图 3.8(a)
axes[1, 2].set_xlabel("fx (cyc/m)")
axes[1, 2].set_title("centered |F| vs analytic (Fig. 3.8a)")
axes[1, 2].legend(["discrete", "analytic"], fontsize=8)

for ax in axes.flat:
    ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_fft_example.png"), dpi=150)
print("已保存 out_fft_example.png")
