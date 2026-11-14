# -*- coding: utf-8 -*-
"""
===== 书 Schmidt PM199 第 1 章 =====
单缝衍射 / 双缝干涉：菲涅耳近场 + 夫琅禾费远场
对应原书式 (1.57)（菲涅耳衍射积分）与式 (1.67)、(1.71)（夫琅禾费衍射）。
重写自用户旧稿 InterferenceDiffraction.m / 双缝干涉和单缝衍射.py：
  - 旧 MATLAB 版 propFresnelTF 漏写了 fft2（H 根本没乘到频谱上），本版已修正；
  - OpenCV 画缝改为解析 rect 组合，去掉 cv2 依赖；
  - 统一采用"胞元中心"采样约定（x 偏移半像素）：40 um 矩形缝恰占 20 个胞元，
    离散谱零点严格落在连续理论位置 m*lambda*z/Dx 上（节点约定会偏差 1+1/M）。
用法：直接运行弹窗显示；`python ch1_slit_diffraction.py <输出目录>` 则保存 PNG。
依赖：numpy + matplotlib（不需要 scipy / cv2）。
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ===== 1. 物理与网格参数 =====
N = 1024            # 网格点数（N x N）
dx = 2e-6           # 采样间隔 2 um（MATLAB 数组从 1 起、Python 从 0 起，坐标式一致）
wvl = 532e-9        # 波长 532 nm（变量名避开 Python 关键字 lambda）
k = 2 * np.pi / wvl

L = N * dx                                    # 网格物理宽度 2.048 mm
x1 = (np.arange(N) - N // 2 + 0.5) * dx       # 胞元中心坐标：-L/2+dx/2 ... L/2-dx/2
X1, Y1 = np.meshgrid(x1, x1)                  # numpy meshgrid 与 MATLAB 同名同义

# 单缝：宽 Dx、高 Dy；双缝：中心距 dx_slit（> Dx，两缝不重叠）
Dx = 40e-6
Dy = 400e-6
dx_slit = 120e-6

U_single = ((np.abs(X1) <= Dx / 2) & (np.abs(Y1) <= Dy / 2)).astype(float)
# 双缝 = 两个 rect 之和：缝中心分别在 x = ±dx_slit/2（相干叠加，直接相加）
U_double = ((np.abs(np.abs(X1) - dx_slit / 2) <= Dx / 2)
            & (np.abs(Y1) <= Dy / 2)).astype(float)

z_near = 5e-3      # 近场观察距离 5 mm（菲涅耳区）
z_far = 80e-3      # 远场观察距离 80 mm（夫琅禾费区）


# ===== 2. 传播算子 =====
def prop_fresnel_tf(U_in, delta, wavelength, z):
    """菲涅耳衍射·传递函数法（TF）：U2 = IFT2{ FT2{U1} * H }，
    H(fx,fy) = exp(i k z) exp(-i pi lambda z (fx^2+fy^2))。
    传递函数法适用条件：pi*lambda*z/(4 dx^2) < pi，即 z < 4 dx^2/lambda（本例 30 mm）。
    频率轴取 fx = (n - N/2)/(N dx)，与胞元中心 x 配对（半像素相移不改变光强）。
    """
    M2, N2 = U_in.shape
    fx = (np.arange(N2) - N2 // 2) / (N2 * delta)
    fy = (np.arange(M2) - M2 // 2) / (M2 * delta)
    FX, FY = np.meshgrid(fx, fy)
    H = np.exp(1j * k * z) * np.exp(-1j * np.pi * wavelength * z * (FX ** 2 + FY ** 2))
    U_ft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(U_in)))   # 中心化谱
    U_out = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(U_ft * H)))
    return U_out


def fraunhofer_fft(U_in, delta, wavelength, z):
    """夫琅禾费衍射（式 1.67）：U2 = e^{ikz} e^{ik(x2^2+y2^2)/(2z)} /(i lambda z) * FT2{U1}，
    观察面坐标由频率缩放得到 x2 = lambda z * fx。返回 (U2, X2, Y2)。
    """
    M2, N2 = U_in.shape
    fx = (np.arange(N2) - N2 // 2) / (N2 * delta)
    fy = (np.arange(M2) - M2 // 2) / (M2 * delta)
    FX, FY = np.meshgrid(fx, fy)
    x2 = wavelength * z * FX          # 远场屏幕真实坐标（米）
    y2 = wavelength * z * FY
    G = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(U_in))) * delta ** 2   # 连续 FT 的离散近似
    U_out = np.exp(1j * k * z) * np.exp(1j * k * (x2 ** 2 + y2 ** 2) / (2 * z)) \
        / (1j * wavelength * z) * G
    return U_out, x2, y2


# ===== 3. 计算 =====
I_single_fres = np.abs(prop_fresnel_tf(U_single, dx, wvl, z_near)) ** 2
I_double_fres = np.abs(prop_fresnel_tf(U_double, dx, wvl, z_near)) ** 2

U_s_far, X2, Y2 = fraunhofer_fft(U_single, dx, wvl, z_far)
U_d_far, _, _ = fraunhofer_fft(U_double, dx, wvl, z_far)
I_single_far = np.abs(U_s_far) ** 2
I_double_far = np.abs(U_d_far) ** 2

x1_um = x1 * 1e6          # 源平面坐标（um）
x2_um = X2[0, :] * 1e6    # 远场横向坐标（um）

# ===== 4. 数值验证（运行后核对打印结果） =====
def padded_far_profile(U_in):
    """4 倍补零后再取远场 |U| 剖面：把极值定位精度提高约 4 倍。
    补零在 fftshift 后的一维序列上进行：正坐标样本（前半段）放头部，
    负坐标样本（后半段）接到尾部，保证采样列在 4N 长度上仍等距连续、
    相位参考点（x≈0）不动，否则频谱会出现虚假调制。
    """
    Np = 4 * N
    h = np.fft.fftshift(U_in[N // 2, :])           # y=0 行
    hp = np.zeros(Np)
    hp[:N // 2] = h[:N // 2]                       # x = +dx/2 ... +L/2-dx/2
    hp[-(N // 2):] = h[N // 2:]                    # x = -L/2+dx/2 ... -dx/2（环形接回尾部）
    F = np.fft.fft(hp) * dx                        # 连续 FT 的离散近似（一维）
    f = np.arange(Np) / (Np * dx)                  # k >= 0 半轴频率
    f[f > 1 / (2 * dx)] -= 1 / dx                  # 折回负频率
    return np.abs(F) ** 2, f * wvl * z_far * 1e6   # (I, x2/um)


print('===== 夫琅禾费衍射数值验证（4 倍补零，y2 = 0 剖面）=====')
I_pad_d, x2_pad = padded_far_profile(U_double)
I_pad_s, _ = padded_far_profile(U_single)

# (1) 左右对称性：在镜像点 -x2 处插值取值再比较。
# 注意：偶数 N 的中心化频率轴是 -N/2...N/2-1 个 bin，本身左右差一个 bin，
# 直接 I 与 I[::-1] 翻转比较会错位一个 bin，必须插值到 -x2。
prof_d = I_double_far[N // 2, :]
mirror = np.interp(-x2_um, x2_um, prof_d)
sym = np.max(np.abs(prof_d - mirror)) / prof_d.max()
print(f'(1) 双缝远场剖面左右对称性  max|I(x2)-I(-x2)|/max(I) = {sym:.2e}')

# (2) 单缝暗纹：sinc 零点（严格）理论位置 x2 = m * lambda*z/Dx
dark_th = np.arange(1, 8) * wvl * z_far / Dx * 1e6
win = (x2_pad > 0) & (x2_pad < dark_th[-1] * 1.1)
xm, pm = x2_pad[win], I_pad_s[win]
dark_num = np.array([xm[i] for i in range(1, len(pm) - 1)
                     if pm[i] < pm[i - 1] and pm[i] < pm[i + 1]])
rel_dark = np.max(np.abs(dark_num[:len(dark_th)] - dark_th) / dark_th)
print(f'(2) 单缝暗纹理论 (um) = {np.round(dark_th, 1)}')
print(f'    数值极小   (um) = {np.round(dark_num[:len(dark_th)], 1)}  最大相对偏差 {rel_dark:.2e}')

# (3) 双缝：暗纹零点严格位于 (n+1/2)*lambda*z/dx_slit 与包络零点 m*lambda*z/Dx；
#     亮纹极大受单缝包络斜坡影响，会从 cos^2 极大 n*lambda*z/dx_slit 向中心偏移
#     （第一亮纹约 341 um 而非 354.7 um），故应与式 (1.71) 解析剖面的极大对比。
cos_zero_th = (np.arange(6) + 0.5) * wvl * z_far / dx_slit * 1e6
env_zero_th = np.arange(1, 7) * wvl * z_far / Dx * 1e6
sub = (x2_pad > 0) & (x2_pad < 6.5 * wvl * z_far / dx_slit * 1e6)
xd, pd = x2_pad[sub], I_pad_d[sub]     # xd/pd 是子集，之后一律用子集索引
# --- 3a 暗纹：取数值剖面在理论零点处的值，与邻近干涉条纹峰值之比 ---
frac = []
sp = wvl * z_far / dx_slit * 1e6          # cos^2 因子的条纹间距 (um)
for z0 in np.concatenate([cos_zero_th[:3], env_zero_th[:2]]):
    i0 = int(np.argmin(np.abs(xd - z0)))
    win = np.abs(xd - z0) < 0.9 * sp      # 邻近条纹峰：零点两侧一个条纹周期内取最大
    frac.append(pd[i0] / pd[win].max())
print(f'(3a) 双缝暗纹深度：理论零点处 I/邻近条纹峰 = {np.round(np.array(frac), 5)}（全部 << 1）')
# --- 3b 亮纹：解析式在细网格上找极大；数值谱直接取 bin 极大 ---
x2_fine = np.arange(0.0, 6.2 * wvl * z_far / dx_slit * 1e6, 0.2)
pa_fine = (np.cos(np.pi * dx_slit * (x2_fine * 1e-6) / (wvl * z_far)) ** 2
           * np.sinc(Dx * (x2_fine * 1e-6) / (wvl * z_far)) ** 2)
bright_ana = np.array([x2_fine[i] for i in range(1, len(pa_fine) - 1)
                       if pa_fine[i] > pa_fine[i - 1] and pa_fine[i] > pa_fine[i + 1]])
bright_num = np.array([xd[i] for i in range(1, len(pd) - 1)
                       if pd[i] > pd[i - 1] and pd[i] > pd[i + 1]])
nb = min(len(bright_num), len(bright_ana))
rel_bright = np.max(np.abs(bright_num[:nb] - bright_ana[:nb]) / bright_ana[:nb])
print(f'(3b) 亮纹 cos^2 极大（朴素）(um) = {np.round(np.arange(1, 7) * sp, 1)}')
print(f'     解析式(1.71)极大     (um) = {np.round(bright_ana[:nb], 1)}  <- 被包络拉向中心')
bin_um = wvl * z_far / (4 * N * dx) * 1e6
print(f'     数值极大             (um) = {np.round(bright_num[:nb], 1)}  （bin 距 {bin_um:.1f} um 量化）')
print(f'     数值 vs 解析极大最大相对偏差 = {rel_bright:.2e}')
miss = [n for n in (3, 6) if np.min(np.abs(bright_num - n * sp)) > 30]
print(f'     缺级检验：n=3,6 处无极大 -> 缺失级次 = {miss}（理论预期 dx_slit/Dx=3 -> [3, 6]）')

# (4) 与式 (1.71) 解析解对比：|U| 沿 y2 = 0，取显示窗口 ±800 um（未补零原网格）
U_ana = 2 * Dx * Dy / (wvl * z_far) * np.cos(np.pi * dx_slit * X2 / (wvl * z_far)) \
    * np.sinc(Dx * X2 / (wvl * z_far)) * np.sinc(Dy * Y2 / (wvl * z_far))
win4 = np.abs(x2_um) <= 800
err = np.max(np.abs(np.abs(U_d_far[N // 2, win4]) - np.abs(U_ana[N // 2, win4]))) \
    / np.abs(U_d_far[N // 2, :]).max()
print(f'(4) FFT 数值解 vs 式(1.71)解析解（±800 um 内）  max||U|-|U_ana||/max|U| = {err:.2e}')

# ===== 5. 绘图 =====
def make_figures():
    fig1, axes = plt.subplots(2, 2, figsize=(11, 9))
    panels = [
        (axes[0, 0], I_single_fres, x1_um, f'单缝·菲涅耳近场 z={z_near*1e3:.0f} mm'),
        (axes[0, 1], I_single_far, x2_um, f'单缝·夫琅禾费远场 z={z_far*1e3:.0f} mm'),
        (axes[1, 0], I_double_fres, x1_um, f'双缝·菲涅耳近场 z={z_near*1e3:.0f} mm'),
        (axes[1, 1], I_double_far, x2_um, f'双缝·夫琅禾费干涉 z={z_far*1e3:.0f} mm'),
    ]
    for ax, I, xc, title in panels:
        im = ax.imshow(I, extent=[xc[0], xc[-1], xc[0], xc[-1]], cmap='gray',
                       vmin=0, vmax=0.02 * I.max(), origin='lower')
        ax.set_title(title)   # vmax 截到 2%：压掉 0 级亮斑，突出条纹
        ax.set_xlabel('x (um)')
        ax.set_ylabel('y (um)')
        ax.set_xlim(-800, 800)
        ax.set_ylim(-800, 800)
        fig1.colorbar(im, ax=ax, shrink=0.85)
    fig1.tight_layout()

    fig2, ax = plt.subplots(figsize=(10, 5))
    norm = I_pad_d.max()
    ax.plot(x2_pad, I_pad_d / norm, 'b-', lw=1.2, label='双缝 FFT 数值解（4 倍补零）')
    env = np.sinc(Dx * (x2_pad * 1e-6) / (wvl * z_far)) ** 2    # 单缝包络 sinc^2
    ax.plot(x2_pad, env, 'g--', lw=1.2, label='单缝包络 sinc$^2$($D_x x_2/\\lambda z$)')
    ax.set_xlim(-800, 800)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel('远场坐标 $x_2$ (um)')
    ax.set_ylabel('归一化光强')
    ax.set_title('双缝夫琅禾费干涉：$\\Delta x = 3D_x$ → ±3 级缺级')
    ax.legend()
    ax.grid(alpha=0.3)
    fig2.tight_layout()
    return fig1, fig2


if __name__ == '__main__':
    if len(sys.argv) > 1:
        outdir = sys.argv[1]
        os.makedirs(outdir, exist_ok=True)
        f1, f2 = make_figures()
        f1.savefig(os.path.join(outdir, 'slits_2x2.png'), dpi=150)
        f2.savefig(os.path.join(outdir, 'slits_profile.png'), dpi=150)
        print('已保存：slits_2x2.png, slits_profile.png ->', outdir)
    else:
        make_figures()
        plt.show()
