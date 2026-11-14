# -*- coding: utf-8 -*-
"""书 7.2.3 节:生成 USAF 1951 分辨率测试板的理想图像(替代书配套 PNG)。

原书 coh_image 脚本用 imread 读 250x250 的 'USAF1951B250.png'(书网站下载,
原地址已失效)。本脚本程序化生成同规格的理想图像 USAF1951_py250.png:
  - 250 x 250 像素、黑底(0)白条(255),映射到像面边长 L = 0.3 mm,
    像素间隔 du = 1.2 um(满足书式 7.18 的采样判据)。
  - USAF 1951 频率标尺:组 g、单元 n 的空间频率 f = f00 * 2^(g + n/6),
    锚定 f00 = f(组0,单元1) = 2e5 cyc/m——使"组-2/单元6"恰为
    1e5 cyc/m(周期 10 um),与书 7.2.3 节"周期约 1e-5 m"的判读一致。
    (原书 PNG 的实际缩放未知,书 7.3.3 节判读组-1/单元3约 1.7e5 cyc/m,
    本板同单元为 1.41e5——同一量级,定性结论全部一致,见笔记正文说明。)
  - 每单元 3 条竖直条,条宽 = 1/(2f),条长 = 5 倍条宽(USAF 规格);
    8 倍超采样后盒式降采样(抗混叠):超出奈奎斯特频率 1/(2du) 的最细
    条纹自然退化为灰色——"理想图像里没有的频率,仿真也造不出来"。

输出:USAF1951_py250.png(250x250 灰度)、out_fig7_5.png(显示图,对应书图 7.5)。

MATLAB/Python 差异:
  - MATLAB 侧 make_usaf_chart.m 生成同规格 PNG;书原脚本 imread('USAF1951B250')
    的读法在 coh_image.m 中改为读本文件生成的 USAF1951_py250。
  - 图像阵列"顶行在前"的存储约定与书一致:读入后 flipud,让图像底部
    对应负 v 坐标(书第 5 行 A=flipud(A))。
"""
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

N_PX = 250                 # 理想图像边长(样本数)
L_IMG = 0.3e-3             # 像面边长 (m)
DU = L_IMG / N_PX          # 像素间隔 = 1.2 um
F00 = 2.0e5                # (组0,单元1) 的空间频率 [cyc/m]
SS = 8                     # 超采样倍数(抗混叠)
GROUPS = (-2, -1, 0, 1, 2) # 组
ELEMENTS = (1, 2, 3, 4, 5, 6)

_HERE = os.path.dirname(os.path.abspath(__file__))
PNG_PATH = os.path.join(_HERE, "USAF1951_py250.png")


def _font(size):
    for name in ("arial.ttf", "Arial.ttf", "segoeui.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def make_chart(path=PNG_PATH):
    """绘制 USAF 式分辨率板并保存为 250x250 灰度 PNG。"""
    size = N_PX * SS
    img = Image.new("L", (size, size), 0)
    d = ImageDraw.Draw(img)
    font_g = _font(13 * SS)   # 组标签
    font_e = _font(9 * SS)    # 单元编号
    x_left = 30 * SS          # 左边距(单元编号)
    y_top = 18 * SS           # 顶边距(组标签)
    col_w = (N_PX - 30) / 5 * SS

    # 组标签(顶行)
    for i, g in enumerate(GROUPS):
        txt = f"{g:+d}" if g > 0 else f"{g}"
        bb = d.textbbox((0, 0), txt, font=font_g)
        tw = bb[2] - bb[0]
        cx = x_left + (i + 0.5) * col_w
        d.text((cx - tw / 2, 2 * SS), txt, fill=255, font=font_g)

    # 逐组逐单元画 3 条竖条
    for i, g in enumerate(GROUPS):
        y = y_top
        x0 = x_left + i * col_w
        for n in ELEMENTS:
            f = F00 * 2.0 ** (g + n / 6.0)     # 该单元空间频率 [cyc/m]
            w = 1.0 / (2 * f) / DU * SS        # 条宽(超采样像素)
            ln = 5 * w                         # 条长 = 5 倍条宽
            if y + ln > size:
                break
            for b in range(3):                 # 3 条,间隔 = 条宽(周期 2w)
                bx = x0 + b * 2 * w
                d.rectangle([bx, y, bx + w, y + ln], fill=255)
            if g == -2:                        # 最左组旁标单元编号(示意)
                d.text((x0 - 24 * SS, y + ln / 2 - 6 * SS),
                       str(n), fill=255, font=font_e)
            y += ln + 4 * SS
    img = img.resize((N_PX, N_PX), Image.BOX)  # 盒式降采样 = 抗混叠
    arr = np.asarray(img, dtype=float) / 255.0
    arr = (arr > 0.5).astype(float)   # 二值化:与原书 PNG 一致(书 7.2.3 节
    #                                # "归一化后只含 0 和 1",sqrt(Ig) 不起作用);
    #                                # 超出奈奎斯特频率的细条纹(灰度<0.5)就此消失,
    #                                # 不给理想图像混入假频率
    Image.fromarray((arr * 255).astype(np.uint8)).save(path)
    return arr
    return np.asarray(img, dtype=float) / 255.0


def ensure_chart(path=PNG_PATH):
    """返回理想图像 Ig(0..1);PNG 不存在时先生成。"""
    if not os.path.exists(path):
        make_chart(path)
    return np.asarray(Image.open(path).convert("L"), dtype=float) / 255.0


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Ig = make_chart()
    # ---- 数值自检:组-2/单元6 的条纹周期应约 10 um ----
    u = (np.arange(N_PX) - N_PX // 2) * DU            # 居中坐标 (m)
    # 组-2 在最左列;跳过顶部组标签行(前 18 px),从第一个单元条处扫一行:
    row = None
    for r in range(18, N_PX):
        seg = Ig[r, 30:74]
        if seg.max() > 0.9 and seg.min() < 0.1:
            row = r
            break
    prof = Ig[row, :]
    # 数过峰(>0.5)的位置(相邻条中心间距 = 空间周期 1/f)
    peaks = []
    for k in range(30, 74):
        if prof[k] > 0.5 and prof[k - 1] <= 0.5:
            peaks.append(u[k])
    if len(peaks) >= 2:
        periods = np.diff(peaks)
        print(f"组-2/单元1 条峰间距     = {np.mean(periods)*1e6:.2f} um "
              f"(理论 1/(2^(-11/6)*2e5) = {1/(F00*2**(-11/6))*1e6:.2f} um)")
    print(f"理想图像灰度范围        = [{Ig.min():.3f}, {Ig.max():.3f}]")
    print(f"图像均值(条形占空比)  = {Ig.mean():.4f}")

    # ---- 显示图(对应书图 7.5) ----
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    ax.imshow(Ig, extent=[u[0], u[-1] + DU, u[0], u[-1] + DU],
              origin="upper", cmap="gray")           # origin=upper:顶行在前
    ax.set_xlabel("u (m)"); ax.set_ylabel("v (m)")
    ax.set_title("USAF1951_py250 ideal image")
    fig.tight_layout()
    fig.savefig(os.path.join(_HERE, "out_fig7_5.png"), dpi=110)
    print("已保存 USAF1951_py250.png / out_fig7_5.png")
