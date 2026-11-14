# ===== 书 1.1 例 1.1 · Table 1.2 Plot_sinc.m 的 Python 版 =====
# 绘制 sinc(x) = sin(pi*x)/(pi*x) 曲线，验证其零点在 x = ±1, ±2, ±3, ...
#
# MATLAB / Python 差异提示：
# 1) MATLAB 的 x=-5.5:0.01:5.5 是闭区间等差序列；numpy.arange 不含终点，
#    这里用 np.linspace 精确复刻同样的 1101 个采样点。
# 2) MATLAB 中 0/0 得到 NaN 并在 plot 里留下一个针孔；numpy 中是 nan + RuntimeWarning。
#    原书代码没处理 x=0，本版在 x=0 处显式补上极限值 1（sinc 的定义值）。

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-5.5, 5.5, 1101)          # 与书中 x=-5.5:0.01:5.5 相同的 1101 点
sinc = np.sin(np.pi * x) / (np.pi * x)    # 原书写法：sinc=sin(pi*x)./(pi*x)
sinc[x == 0] = 1.0                        # 修补 x=0 处的 0/0（原书未处理，见笔记"坑与陷阱"）

plt.plot(x, sinc)
plt.axis([-5.5, 5.5, -0.3, 1.1])
plt.grid(True)
plt.xlabel('x')
plt.ylabel('sinc (x)')
plt.title('sinc function')
plt.show()
