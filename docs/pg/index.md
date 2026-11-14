---
icon: material/gamepad-variant
---

# Python极客项目编程

Mahesh Venkitachalam《Python Playground》全书 **14 章**精读笔记，五个板块：热身运动（iTunes 解析、万花尺）→ 模拟生命（生命游戏、Karplus-Strong 弦音、类鸟群）→ 图片之乐（ASCII 艺术、照片马赛克、三维立体画）→ 走进三维（OpenGL、粒子烟花、体渲染）→ 玩转硬件（Arduino 光敏、激光音乐秀、树莓派气象站）。

!!! note "运行环境与代码策略"
    conda 环境 `machinelearn`（Python 3.10 + NumPy + Matplotlib + Pillow）。原书基于 Python 3.3，笔记代码已全部核对为 Python 3 可运行写法；绘图 / 模拟 / 图像类项目（第 1、3–8 章）在本机实跑出结果图（标注"代码实跑"），`turtle` 窗口类项目用 Matplotlib 复现数学曲线，OpenGL 与硬件类项目（第 9–14 章）忠实呈现代码与电路原理，概念示意图截自原书。每章的分节代码在 [nb/](nb/) 目录（`chNN.ipynb`，含 Arduino sketch cell）；第 1、3–8 章带本机固定种子实跑输出，第 9–14 章为完整实现。

<div class="grid cards" markdown>

-   **第1章 · 解析 iTunes 播放列表**

    ---

    plistlib 读 XML 播放列表、集合找重复音轨与共同音轨、matplotlib 直方图与散点图。

    [:octicons-arrow-right-24: 开始](ch01.md)

-   **第2章 · 万花尺**

    ---

    参数方程与海龟作图：内旋轮线方程、GCD 定周期、Spiro 类与动画。

    [:octicons-arrow-right-24: 开始](ch02.md)

-   **第3章 · Conway 生命游戏**

    ---

    元胞自动机四规则、环形边界条件、numpy 网格与 matplotlib 动画、滑翔机图案。

    [:octicons-arrow-right-24: 开始](ch03.md)

-   **第4章 · 用 Karplus-Strong 算法产生音乐泛音**

    ---

    环形缓冲区模拟拨弦、低通滤波与衰减、deque 实现、WAV 文件与五声音阶。

    [:octicons-arrow-right-24: 开始](ch04.md)

-   **第5章 · 类鸟群：仿真鸟群**

    ---

    Reynolds 三规则（分离/列队/内聚）、向量化 numpy 计算提速 200 倍、驱散交互。

    [:octicons-arrow-right-24: 开始](ch05.md)

-   **第6章 · ASCII 文本图形**

    ---

    灰度亮度映射字符梯度、平均亮度计算、横纵比校正、10 级与 70 级灰度。

    [:octicons-arrow-right-24: 开始](ch06.md)

-   **第7章 · 照片马赛克**

    ---

    目标图像网格分割、平均 RGB 匹配、距离比较选块、PIL 拼接成品图。

    [:octicons-arrow-right-24: 开始](ch07.md)

-   **第8章 · 三维立体画**

    ---

    线性间距与深度知觉、深度图、像素水平偏移编码深度、"墙眼"观看原理。

    [:octicons-arrow-right-24: 开始](ch08.md)

-   **第9章 · 理解 OpenGL**

    ---

    可编程图形管线、GLSL 顶点/片段着色器、VAO/VBO、纹理贴图与三维变换。

    [:octicons-arrow-right-24: 开始](ch09.md)

-   **第10章 · 粒子系统**

    ---

    抛物线运动建模、GPU 着色器算粒子轨迹、alpha 混合与公告板、火花喷泉动画。

    [:octicons-arrow-right-24: 开始](ch10.md)

-   **第11章 · 体渲染**

    ---

    体光线投射算法、彩色立方体生成光线、FBO 离屏渲染、三维纹理与二维切片。

    [:octicons-arrow-right-24: 开始](ch11.md)

-   **第12章 · Arduino 简介**

    ---

    光敏电阻分压电路、模拟引脚读取、串口通信、pySerial + matplotlib 实时滚动图。

    [:octicons-arrow-right-24: 开始](ch12.md)

-   **第13章 · 激光音乐秀**

    ---

    旋转镜反射产生图案、FFT 频段分析、PWM 电机控制、音频驱动激光投影。

    [:octicons-arrow-right-24: 开始](ch13.md)

-   **第14章 · 基于树莓派的天气监控器**

    ---

    DHT11 传感器、Bottle Web 框架、flot 浏览器绘图、AJAX 轮询与 GPIO 控制 LED。

    [:octicons-arrow-right-24: 开始](ch14.md)

</div>
