# -*- coding: utf-8 -*-
"""书 7.2.4 节粗糙物体的相干成像(书图 7.10 的 speckle)+ 习题 7.5。

7.2.4:ug = sqrt(Ig) .* exp(j*2*pi*rand(M)) —— 每个样本一个 [0,2pi) 独立
均匀相位,代表相对光波长"极粗糙"的表面;相干成像后得到 speckle(相干斑)。
习题 7.5(a):把 2pi 乘子调小(相位起伏幅度下降),像里的斑纹对比度下降,
图案重新可辨——相位差越小,各点波前越接近"同相",干涉起伏越弱。
习题 7.5(b):改孔径(f/5 vs f/20):PSF 变宽 -> 斑纹颗粒变大
(颗粒线度约 1.22*lam*(f/#)),与衍射限分辨本领一致。
输出:out_fig7_10.png(书图 7.10 复现)、out_fig7_ex75.png((a)(b)两栏)
+ 命令行数值检查(斑纹对比度、颗粒自相关宽度)。

MATLAB/Python 差异:
  - MATLAB rand(M) 与 numpy default_rng().random((M,M)) 都是 [0,1) 均匀
    随机数,但种子/序列不同——复现的 speckle 图样与书图"统计一致、
    逐像素不同",这是随机相位仿真的固有性质。
  - exp(j*2*pi*rand) 在 numpy 写作 np.exp(1j*2*np.pi*rng.random((M,M)))。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
sys.path.insert(0, os.path.join(_here, "..", "common"))
from base import coords                                  # noqa: E402
from ch7_imaging import freq_coords, H_circ, coh_image   # noqa: E402
from ch7_usaf_chart import ensure_chart, L_IMG, N_PX     # noqa: E402

rng = np.random.default_rng(7)       # 固定种子保证笔记数值可复现

Ig = ensure_chart()
Ig = np.flipud(Ig)
Ig = Ig / Ig.max()

L, M = L_IMG, N_PX
du = L / M
u, _ = coords(L, M)
lam = 0.5e-6
wxp, zxp = 6.25e-3, 125e-3
f0 = wxp / (lam * zxp)
fu = freq_coords(L, M)
H = H_circ(fu, f0)


def coh_rough(phase_scale, fn=10.0):
    """粗糙表面相干成像:phase_scale 为随机相位的 [0, phase_scale] 幅度。

    phase_scale = 2*pi 即书 7.2.4 的 exp(j*2*pi*rand);调小即习题 7.5(a)。
    """
    ug = np.sqrt(Ig) * np.exp(1j * phase_scale * rng.random((M, M)))
    Hk = H_circ(freq_coords(L, M), 1.0 / (lam * fn) / 2)
    ui, _ = coh_image(ug, Hk)
    return np.abs(ui) ** 2


def speckle_contrast(Ii, mask=None):
    """斑纹对比度 C = std(Ii)/mean(Ii)。

    充分发展斑纹(圆复高斯场)的辐照度服从指数分布,理论 C = 1。
    mask 给出"被照明"的像素区域(按理想图像膨胀得到),排除全黑背景。"""
    roi = Ii if mask is None else Ii[mask]
    return np.std(roi) / np.mean(roi)


def grain_width(Ii):
    """斑纹颗粒线度:Ii 起伏的自相关第一过零距离(沿 u)。

    只对"均匀照明 + 全场粗糙"的纯斑纹图像用(没有物体结构干扰):
    理论颗粒线度与 PSF 同量级,约 1.22*lam*(f/#)。"""
    a = Ii.mean() * np.ones_like(Ii)
    fl = np.fft.fft2(Ii - a)
    ac = np.fft.ifftshift(np.real(np.fft.ifft2(np.abs(fl) ** 2)))
    ac /= ac[M // 2, M // 2]
    line = ac[M // 2, M // 2:]                     # 从相关峰向 +u
    idx = int(np.argmax(line[1:] < 0.02)) + 1      # 主瓣外第一处跌到 ~0
    return idx * du


def coh_rough_uniform(phase_scale, fn):
    """均匀振幅 + 随机相位的"全粗糙"物体(纯斑纹,用于定量)。"""
    ug = np.exp(1j * phase_scale * rng.random((M, M)))
    Hk = H_circ(freq_coords(L, M), 1.0 / (lam * fn) / 2)
    ui, _ = coh_image(ug, Hk)
    return np.abs(ui) ** 2


# ---------- 书 7.2.4(书图 7.10):2*pi 全粗糙,f/10 ----------
Ii_full = coh_rough(2 * np.pi, 10.0)

fig, ax = plt.subplots(figsize=(4.6, 4.2))
ax.imshow(Ii_full, extent=[u[0], u[-1] + du, u[0], u[-1] + du],
          origin="lower", cmap="gray")
ax.set_xlabel("u (m)"); ax.set_ylabel("v (m)")
ax.set_title("coherent image of a rough object (speckle)")
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_10.png"), dpi=110)
plt.close(fig)

# ---------- 习题 7.5(a):粗糙度递减; (b):孔径变化 ----------
# 上排:测试板 + 相位尺度递减(书题意:图像上斑纹变弱、图案重新可辨);
# 下排:均匀粗糙物体的纯斑纹(f/5 -> f/20,颗粒随 f/# 变大)。
# 定量对比度用均匀物体测(测试板条形结构会污染 std/mean,只作视觉)。
fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.6))
scales = [2 * np.pi, np.pi, 0.4 * np.pi]
for k, s in enumerate(scales):
    Ii = coh_rough(s, 10.0)
    axs[0, k].imshow(np.sqrt(Ii), extent=[u[0], u[-1] + du, u[0], u[-1] + du],
                     origin="lower", cmap="gray")
    axs[0, k].set_title(f"chart, phase scale = {s/np.pi:.1f}pi")
    axs[0, k].set_xticks([]); axs[0, k].set_yticks([])
for k, fn in enumerate((5.0, 10.0, 20.0)):
    Ii_u = coh_rough_uniform(2 * np.pi, fn)
    axs[1, k].imshow(np.sqrt(Ii_u), extent=[u[0], u[-1] + du, u[0], u[-1] + du],
                     origin="lower", cmap="gray")
    axs[1, k].set_title(f"uniform rough, f/{fn:g}, grain = "
                        f"{grain_width(Ii_u)*1e6:.1f} um")
    axs[1, k].set_xticks([]); axs[1, k].set_yticks([])
fig.suptitle("Ex 7.5: (top) roughness down -> speckle contrast down; "
             "(bottom) aperture down -> speckle grains bigger", y=0.995)
fig.tight_layout(); fig.savefig(os.path.join(_here, "out_fig7_ex75.png"), dpi=105)
plt.close(fig)

print("== 书 7.2.4 + 习题 7.5 数值检查 ==")
Ii_u = coh_rough_uniform(2 * np.pi, 10.0)
print(f"均匀粗糙 2*pi(f/10)      : C = {speckle_contrast(Ii_u):.3f}"
      f"(充分发展斑纹理论 C = 1)")
for s in scales:
    Ii_u = coh_rough_uniform(s, 10.0)
    print(f"相位尺度 {s/np.pi:4.1f}*pi(f/10): C = {speckle_contrast(Ii_u):.3f}")
for fn in (5.0, 10.0, 20.0):
    Ii_u = coh_rough_uniform(2 * np.pi, fn)
    print(f"f/{fn:g} 均匀粗糙 2*pi      : 颗粒自相关过零 = {grain_width(Ii_u)*1e6:5.1f} um"
          f"(1.22*lam*(f/#) = {1.22*lam*fn*1e6:.1f} um), C = "
          f"{speckle_contrast(Ii_u):.2f}")
print("已保存 out_fig7_10.png / out_fig7_ex75.png")
