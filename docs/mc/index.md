---
icon: material/dice-multiple
---

# 随机模拟方法及应用

肖柳青、周石鹏《随机模拟方法与应用》（北京大学出版社，2014）全书 **11 章**精读笔记：从蒲丰投针和蒙提霍尔问题入手，经概率论复习、随机数生成与抽样，到 MCMC、排队模拟、蒙特卡罗优化、随机游走、蒙特卡罗积分与复杂系统模拟。

!!! note "知识来源与章节顺序"
    章节顺序完全跟随中文教材（11 章，线性推进，后面内容不提前）——方便对照教材逐章学习。知识构建主要依据两本英文教材：**Dunn & Shultis《Exploring Monte Carlo Methods》(2012)** 与 **Shonkwiler & Mendivil《Explorations in Monte Carlo Methods》(2e)**，理论推导以这两本为准并适当加墨，顺带复习概率论与数理统计。中文教材负责提纲和例子，英文教材负责原理和推导。

!!! note "运行环境约定"
    Python 3.12（Anaconda base）+ NumPy + Matplotlib + SciPy。教材用 Matlab 的地方，笔记一律换成 **Python** 实现，差异点在正文标注。每章代码实现在本板块 `nb/` 目录下的 Jupyter Notebook（页面内代码块与 notebook 一致，notebook 可下载复跑）；结果插图均来自本机实跑（标注"代码实跑"）。

<div class="grid cards" markdown>

-   **第1章 · 初识随机模拟方法**

    ---

    蒲丰投针、蒙提霍尔、优惠券收集、投点求 π 与定积分——五个经典例子跑通"模拟是什么"。

    [:octicons-arrow-right-24: 开始](ch01.md)

-   **第2章 · 懂点概率论**

    ---

    概率公理、条件概率与 Bayes 公式、随机变量与分布、期望方差、变换、大数定律与中心极限定理（推导不跳步）。

    [:octicons-arrow-right-24: 开始](ch02.md)

-   **第3章 · 善用身边的数学秘书**

    ---

    教材讲 Matlab，这里换成 Python/NumPy/Matplotlib/Jupyter 全套对照，另附随机模拟常用函数速查。

    [:octicons-arrow-right-24: 开始](ch03.md)

-   **第4章 · 让电脑玩掷骰子：生成随机数**

    ---

    平方取中、线性同余生成器、周期与谱检验、组合生成器、Mersenne Twister 与随机数统计检验。

    [:octicons-arrow-right-24: 开始](ch04.md)

-   **第5章 · 掷骰子的进阶：特殊分布抽样**

    ---

    逆变换法、舍选法、组合法、alias 方法、Box-Muller、指数与泊松过程、Gamma/Beta 抽样。

    [:octicons-arrow-right-24: 开始](ch05.md)

-   **第6章 · 神奇的马尔科夫链蒙特卡罗方法**

    ---

    马尔可夫链、平稳分布与细致平衡、Metropolis 与 Metropolis-Hastings、Gibbs 抽样、贝叶斯推断应用。

    [:octicons-arrow-right-24: 开始](ch06.md)

-   **第7章 · 仿真随机服务系统**

    ---

    离散事件模拟（事件调度法）、M/M/1 排队、Little 定律的模拟验证。

    [:octicons-arrow-right-24: 开始](ch07.md)

-   **第8章 · 不落俗套：蒙特卡罗优化方法**

    ---

    模拟退火与冷却调度、TSP 与路径反转、遗传算法。

    [:octicons-arrow-right-24: 开始](ch08.md)

-   **第9章 · 模拟醉汉行走：随机游走模拟**

    ---

    1d 随机游走、扩散方程、布朗运动、赌徒破产、几何布朗运动与期权定价。

    [:octicons-arrow-right-24: 开始](ch09.md)

-   **第10章 · 化腐朽为神奇：蒙特卡罗积分法**

    ---

    hit-or-miss 与样本均值估计、方差缩减全套：重要性抽样、分层、系统抽样、对偶变量、控制变量、Russian roulette。

    [:octicons-arrow-right-24: 开始](ch10.md)

-   **第11章 · 随机模拟无极限：探索复杂性**

    ---

    元胞自动机交通流、渗流、Ising 相变回顾、复杂网络上的模拟。

    [:octicons-arrow-right-24: 开始](ch11.md)

</div>
