# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 1 章（衔接第 2 章主题）=====
离散傅里叶变换的三步演化：空间采样 -> 有限截断 -> 频域采样。
演示 g(x)=exp(-a|x|) 从连续傅里叶变换到 DFT 的每一步在两个域里各损失/改变了什么。
重写自用户旧稿 Alleviating_Effects_of_Discretization.m / 采样结果.py，参数一致：
a=10, delta=0.0375, L=0.6, N=16。
用法：直接运行弹窗显示；`python ch1_dft_three_steps.py <输出目录>` 则保存 PNG。
依赖：numpy + matplotlib。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ===== 1. 参数 =====
a = 10.0                    # 指数衰减系数
delta = 0.0375              # 空间采样间隔（MATLAB: delta；Python 避开关键字取名 delta）
L = 0.6                     # 截断窗口宽度
N = int(round(L / delta))   # 采样点数 N = 16

# ===== 2. 连续参考 =====
N_cont = 2000
x_cont = np.linspace(-1.5, 1.5, N_cont)
g_cont = np.exp(-a * np.abs(x_cont))
f_cont = np.linspace(-25.0, 25.0, N_cont)
G_cont = 2 * a / (a ** 2 + (2 * np.pi * f_cont) ** 2)   # FT{exp(-a|x|)} 解析式

# ===== 3. 三步演化 =====
# --- 步骤 1：仅采样（空间域乘 comb，频域卷积 comb -> 频谱以 1/delta 为周期复制）---
x_samp = np.arange(-1.5, 1.5 + delta / 2, delta)
g_samp = np.exp(-a * np.abs(x_samp))
G_repl = np.zeros_like(f_cont)
for p in range(-2, 3):                       # 叠 5 个相邻副本示意混叠
    G_repl += 2 * a / (a ** 2 + (2 * np.pi * (f_cont - p / delta)) ** 2)

# --- 步骤 2：有限截断（空间域乘 rect(x/L)，频域卷积 L*sinc(L f) -> 波纹）---
x_trunc = (np.arange(N) - N // 2) * delta
g_trunc = np.exp(-a * np.abs(x_trunc))
M = 4096                                     # 高倍补零，模拟"连续的带波纹频谱"
g_pad = np.zeros(M)
g_pad[M // 2 - N // 2: M // 2 + N // 2] = g_trunc
G_rip = np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(g_pad)))) * delta
f_dense = (np.arange(M) - M // 2) / (M * delta)

# --- 步骤 3：频域采样（最终 DFT：频域乘 comb(L f)，空间域周期延拓）---
df = 1.0 / (N * delta)
f_dft = (np.arange(N) - N // 2) * df
G_dft = np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(g_trunc)))) * delta

# ===== 4. 数值验证 =====
print('===== DFT 三步演化数值验证 =====')
# (1) 混叠量：步骤 1 的复制谱在 f=1/(2 delta) 处相对连续谱的抬升
f_nyq = 1 / (2 * delta)
ali = (G_repl[np.argmin(np.abs(f_cont - f_nyq))] / G_cont[np.argmin(np.abs(f_cont - f_nyq))])
print(f'(1) 奈奎斯特频率 {f_nyq:.2f} 处 复制谱/连续谱 = {ali:.3f}（采样间隙造成的频谱抬升）')

# (2) 最终 DFT vs 连续谱：在 DFT 频点 f_m = m/(N delta) 上比较
G_dft_bin = np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(g_trunc)))) * delta
G_ana_bin = 2 * a / (a ** 2 + (2 * np.pi * f_dft) ** 2)
relerr = np.max(np.abs(G_dft_bin - G_ana_bin)) / G_ana_bin.max()
print(f'(2) DFT 谱线 vs 连续解析谱（同频点）最大偏差 = {relerr:.2e}'
      f'（来源：截断丢尾 + 邻副本混叠，无法归零）')

# (3) Parseval 能量守恒（DFT 对自身的精确性质，检查实现无误）
lhs = np.sum(g_trunc ** 2) * delta
rhs = np.sum((np.fft.fft(np.fft.fftshift(g_trunc)) * delta) ** 2).real * df
print(f'(3) Parseval：sum|g|^2*delta = {lhs:.6e}，sum|G|^2*df = {rhs:.6e}，'
      f'相对偏差 {abs(lhs-rhs)/lhs:.2e}')

# ===== 5. 绘图 =====
def make_figures():
    fig, axes = plt.subplots(3, 2, figsize=(12, 10))

    axes[0, 0].plot(x_cont, g_cont, color='0.75', lw=1.5)
    axes[0, 0].stem(x_samp, g_samp, linefmt='r-', markerfmt=' ', basefmt='r-')
    axes[0, 0].set_xlim(-1.2, 1.2); axes[0, 0].set_ylim(0, 1.1)
    axes[0, 0].set_title('步骤 1 空间采样：g(x)·comb(x/δ)（无限采样列）')
    axes[0, 0].set_ylabel('g(x)')

    axes[0, 1].plot(f_cont, G_cont, color='0.75', lw=1.5)
    axes[0, 1].plot(f_cont, G_repl, 'r-', lw=2)
    axes[0, 1].set_xlim(-20, 20); axes[0, 1].set_ylim(0, 0.25)
    axes[0, 1].set_title('步骤 1 频域：卷积 comb(δf) → 频谱周期复制（混叠）')
    axes[0, 1].set_ylabel('G(f)')

    out = np.abs(x_samp) >= L / 2
    axes[1, 0].vlines(x_samp[out], 0, g_samp[out], color='0.75', lw=1.2)  # 被丢弃的采样点
    axes[1, 0].stem(x_trunc, g_trunc, linefmt='b-', markerfmt='bo', basefmt='b-')
    axes[1, 0].set_xlim(-1.2, 1.2); axes[1, 0].set_ylim(0, 1.1)
    axes[1, 0].set_title('步骤 2 有限截断：·rect(x/L)（窗口外采样点丢弃）')
    axes[1, 0].set_ylabel('g(x)')

    axes[1, 1].plot(f_cont, G_repl, color='0.75', lw=1.5)
    axes[1, 1].plot(f_dense, G_rip, 'b-', lw=2)
    axes[1, 1].set_xlim(-20, 20); axes[1, 1].set_ylim(0, 0.25)
    axes[1, 1].set_title('步骤 2 频域：卷积 L·sinc(Lf) → 吉布斯波纹')
    axes[1, 1].set_ylabel('G(f)')

    for off in (-L, 0, L):
        axes[2, 0].plot(x_cont + off, g_cont, '--', color='0.75', lw=1)
    axes[2, 0].stem(x_trunc, g_trunc, linefmt='g-', markerfmt='go', basefmt='g-')
    axes[2, 0].set_xlim(-1.2, 1.2); axes[2, 0].set_ylim(0, 1.1)
    axes[2, 0].set_title('步骤 3 频域采样：空间域卷积 comb(x/L) → 周期延拓')
    axes[2, 0].set_ylabel(r'$\tilde{g}(x)$')

    axes[2, 1].plot(f_dense, G_rip, color='0.75', lw=1.5)
    axes[2, 1].stem(f_dft, G_dft, linefmt='g-', markerfmt='go', basefmt='g-')
    axes[2, 1].set_xlim(-20, 20); axes[2, 1].set_ylim(0, 0.25)
    axes[2, 1].set_title('步骤 3 频域：×comb(Lf) → 最终 N 点 DFT 谱线')
    axes[2, 1].set_ylabel(r'$\tilde{G}(f)$')

    for ax in axes.ravel():
        ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig


if __name__ == '__main__':
    if len(sys.argv) > 1:
        outdir = sys.argv[1]
        os.makedirs(outdir, exist_ok=True)
        fig = make_figures()
        fig.savefig(os.path.join(outdir, 'dft_three_steps.png'), dpi=150)
        print('已保存：dft_three_steps.png ->', outdir)
    else:
        make_figures()
        plt.show()
