# ===== 书 3.6 节 + 习题 3.2：两个高斯函数的 FFT 卷积与解析对比 =====
# 参数（书 3.6 节）：wa = 0.3 m、wb = 0.2 m——注释原文 "[exp(-pi) radius]"，
# 意思是 fa(wa) = exp(-pi) ≈ 0.043 的半径（不是 1/e 半径！）；L = 2、M = 200。
# 支撑判据（式 2.27）：Da ≈ 0.9、Db ≈ 0.6，和 1.5 < L = 2 -> 卷绕可忽略。
# 书的技巧：fa、fb 都居中排列，FFT 前【不移位】——两个函数携带同一线性
# 相位，乘积相位对应"整周期平移"，模 L 不可见；但逆变换结果停在"移位
# 排列"（峰值在 index 0），所以要 fftshift 居中。只移位其中一个函数，
# 卷积结果就会整体错位 L/2（本脚本有反例演示）。
# MATLAB/Python 差异：MATLAB 要 .* 与 .^ 表示逐元素，忘点号不报错但结果
# 离谱；numpy 的 * 与 ** 天然逐元素，矩阵乘法反而要另写 @ / np.matmul
# ——两边的"默认值"方向恰好相反。
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from base import coords   # 公共积木：式 (2.21)/(2.22) 的网格坐标（ch2 起）

wa, wb = 0.3, 0.2         # 高斯宽度 [exp(-pi) 半径] (m)
L, M = 2.0, 200           # 边长 (m)、样本数
dx = L / M                # 采样间隔 (m)
x, _ = coords(L, M)       # x 坐标（式 2.21）
fa = np.exp(-np.pi * x ** 2 / wa ** 2)   # 高斯 a
fb = np.exp(-np.pi * x ** 2 / wb ** 2)   # 高斯 b

# ---------- FFT 卷积（书 3.6 节，输入不移位） ----------
Fa = np.fft.fft(fa)                       # 变换 fa
Fb = np.fft.fft(fb)                       # 变换 fb
F0 = Fa * Fb                              # 逐点相乘（numpy 的 * 就是 MATLAB 的 .*）
f0 = np.fft.ifft(F0) * dx                 # 逆变换；卷积定标只乘一次 dx（不是 dx^2）
f = np.fft.fftshift(f0.real)              # 居中（书图 3.9b；虚部是机器精度残渣）

# ---------- 解析卷积（习题 3.2） ----------
# exp(-pi*x^2/w^2) <-> w*exp(-pi*w^2*fx^2)（表 1.3 高斯对）。频域相乘：
# wa*wb*exp(-pi*(wa^2+wb^2)*fx^2)，逆变换回去：
#   c(x) = wa*wb/sqrt(wa^2+wb^2) * exp(-pi*x^2/(wa^2+wb^2))
c_an = wa * wb / np.sqrt(wa ** 2 + wb ** 2) * np.exp(-np.pi * x ** 2 / (wa ** 2 + wb ** 2))

print(f"解析卷积峰值 = {c_an.max():.6f}；FFT 卷积峰值 = {f.max():.6f}")
print(f"FFT 卷积与解析最大差异 = {np.max(np.abs(f - c_an)):.2e}"
      "（高斯衰减极快：离散求和与连续积分之差达机器精度）")

# ---------- 反例：只移位一个函数 -> 卷积整体错位 ----------
fa_s = np.fft.ifftshift(fa)               # 只把 fa 挪到移位排列
f_wrong = np.fft.fftshift((np.fft.ifft(np.fft.fft(fa_s) * np.fft.fft(fb)) * dx).real)
iw = int(np.argmax(f_wrong))
print(f"只移位 fa：卷积峰出现在 x = {x[iw]:.2f} m（应为 0 -> 整体错位 L/2）")

# ---------- 画图：复现书图 3.9 + 解析对比 ----------
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
axes[0].plot(x, fa, label=r"$f_a$ ($w_a$=0.3)")
axes[0].plot(x, fb, "--", label=r"$f_b$ ($w_b$=0.2)")
axes[0].set_xlabel("x (m)"); axes[0].set_title("functions (Fig. 3.9a)")
axes[0].legend(fontsize=8)
axes[1].plot(x, f, label="FFT")
axes[1].plot(x, c_an, ":", label="analytic")
axes[1].set_xlabel("x (m)")
axes[1].set_title("convolution: FFT vs analytic (Fig. 3.9b + Ex. 3.2)")
axes[1].legend(fontsize=8)
for ax in axes:
    ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "out_conv_example.png"), dpi=150)
print("已保存 out_conv_example.png")
