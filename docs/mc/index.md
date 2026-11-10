---
icon: material/dice-multiple
---

# 随机模拟方法及应用

《随机模拟方法及应用》**12 章**精读笔记：从蒲丰投针入手，经概率论基础（复习）、伪随机数与抽样、方差缩减，到 MCMC、贝叶斯推断、随机游走、蒙特卡罗优化、逆问题、线性算子方程与粒子输运模拟。

!!! note "知识来源与章节顺序"
    知识构建主要依赖两本英文教材：**Dunn & Shultis《Exploring Monte Carlo Methods》(2012)** 与 **Shonkwiler & Mendivil《Explorations in Monte Carlo Methods》(2e, 2024)**。章节顺序以 **Dunn & Shultis 章序为主干**，Shonkwiler & Mendivil 的对应内容按主题并入相应章节（每章页首标注对照章节），线性推进、后面内容不提前——方便对照英文教材逐章学习。理论推导以这两本为准并多费笔墨，顺带复习概率论与数理统计（ch02/ch04/ch07 是复习重点章）。

!!! note "运行环境约定"
    Python 3.12（Anaconda base）+ NumPy + SciPy + Matplotlib。每章代码实现在本板块 `nb/` 目录下的 Jupyter Notebook（页面内代码块与 notebook 一致，notebook 可下载复跑）；固定随机种子，结果插图均来自本机实跑（标注"代码实跑"）。

<div class="grid cards" markdown>

-   **第1章 · 初识蒙特卡罗方法**

    ---

    蒲丰投针、频率估概率、求积与模拟两种视角——蒙特卡罗是什么、从哪来、能干什么。[D&S 1 ∪ S&M 1.1/1.4]

    [:octicons-arrow-right-24: 开始](ch01.md)

-   **第2章 · 概率论基础：两大理论支柱**

    ---

    PDF/CDF、期望方差、条件概率与 Bayes、联合分布、大数定律（含弱 LLN 证明）与中心极限定理——概率论复习章。[D&S 2 ∪ S&M 1.2 ∪ 附录B/C]

    [:octicons-arrow-right-24: 开始](ch02.md)

-   **第3章 · 伪随机数**

    ---

    线性同余生成器、满周期判据、谱结构与超平面、组合生成器、随机数统计检验。[D&S 3 ∪ S&M 1.3]

    [:octicons-arrow-right-24: 开始](ch03.md)

-   **第4章 · 抽样、记分与精度**

    ---

    逆变换法、舍选法、组合法、alias 方法，正态/Gamma/Beta 抽样；样本均值与方差、标准误差与置信区间——数理统计复习章。[D&S 4 ∪ S&M 2]

    [:octicons-arrow-right-24: 开始](ch04.md)

-   **第5章 · 方差缩减**

    ---

    重要性抽样、系统抽样、分层抽样、对偶变量、相关抽样、Russian roulette 与 splitting、拟随机数。[D&S 5 ∪ S&M 6.3.1]

    [:octicons-arrow-right-24: 开始](ch05.md)

-   **第6章 · 马尔科夫链蒙特卡罗**

    ---

    遍历马尔科夫链、平稳分布与细致平衡、Metropolis 与 Metropolis-Hastings、Gibbs 抽样、Ising 模型。[D&S 6.1 ∪ S&M 3]

    [:octicons-arrow-right-24: 开始](ch06.md)

-   **第7章 · 贝叶斯推断与 Gibbs 抽样**

    ---

    先验、似然与后验，共轭分布、MAP 与后验预测分布，MCMC 算后验、二维 Gibbs 实例。[D&S 6.2–6.4 ∪ S&M 6.1–6.2]

    [:octicons-arrow-right-24: 开始](ch07.md)

-   **第8章 · 随机游走**

    ---

    1d 随机游走与扩散方程、布朗运动、赌徒破产、几何布朗运动与期权定价、Kelly 准则、电阻网络与 kMC。[S&M 4]

    [:octicons-arrow-right-24: 开始](ch08.md)

-   **第9章 · 蒙特卡罗优化：模拟退火与遗传算法**

    ---

    Metropolis 接受准则与冷却调度、TSP 部分路径反转、遗传算法的编码/交叉/变异。[S&M 5]

    [:octicons-arrow-right-24: 开始](ch09.md)

-   **第10章 · 逆蒙特卡罗**

    ---

    逆问题的三种表述、迭代 IMC、符号蒙特卡罗（SMC）与不确定度量化。[D&S 7]

    [:octicons-arrow-right-24: 开始](ch10.md)

-   **第11章 · 线性算子方程**

    ---

    随机游走解线性方程组与伴随方程、Fredholm 积分方程、Laplace/Poisson 方程的连续蒙特卡罗。[D&S 8]

    [:octicons-arrow-right-24: 开始](ch11.md)

</div>
