---
icon: material/dice-multiple
---

# 随机模拟方法及应用

《随机模拟方法及应用》**12 章**精读笔记：从蒲丰投针入手，经概率论基础（复习）、伪随机数与抽样、方差缩减，到 MCMC、贝叶斯推断、随机游走、蒙特卡罗优化、逆问题、线性算子方程与粒子输运模拟。

!!! note "知识来源与章节顺序"
    以 **Shonkwiler & Mendivil《Explorations in Monte Carlo Methods》(2e, 2024) 章序为主干**，综合全面构建；**Dunn & Shultis《Exploring Monte Carlo Methods》(2012)** 作补充（其独有章节排在其后）。每章页首标注对照章节，线性推进、后面内容不提前——方便对照英文教材逐章学习。理论推导不跳步，顺带复习概率论与数理统计（ch02/ch04/ch08 是复习重点章）。

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

-   **第3章 · 随机数生成**

    ---

    对生成器五要求、平方取中、线性同余与满周期定理、超平面与谱检验、minstd 与组合生成器、卡方与 KS 检验。[S&M 1.3/附录A ∪ D&S 3]

    [:octicons-arrow-right-24: 开始](ch03.md)

-   **第4章 · 概率分布与抽样**

    ---

    CDF 逆变换（离散/连续）、Walker 别名法、二项与泊松、指数与离散事件模拟、正态三招、对数正态与 Gibrat、舍选与 Beta、复合分布与 Gamma、联合分布抽样——数理统计复习章。[S&M 2 ∪ D&S 4]

    [:octicons-arrow-right-24: 开始](ch04.md)

-   **第5章 · 马尔科夫链蒙特卡罗**

    ---

    离散马尔科夫链与平稳分布、Metropolis 与细致平衡、遍历定理、统计力学与 Boltzmann 分布、Ising 模型、Metropolis–Hastings、计数与应用。[S&M 3 ∪ D&S 6.1]

    [:octicons-arrow-right-24: 开始](ch05.md)

-   **第6章 · 随机游走**

    ---

    1d 随机游走与常返性、扩散方程、布朗运动与维纳过程、期权定价、赌徒破产、Kelly 准则、电阻网络与 kMC。[S&M 4]

    [:octicons-arrow-right-24: 开始](ch06.md)

-   **第7章 · 蒙特卡罗优化：模拟退火与遗传算法**

    ---

    Metropolis 接受与冷却调度（Hajek 定理）、TSP 部分路径反转、遗传算法编码/交叉/变异、并行加速。[S&M 5]

    [:octicons-arrow-right-24: 开始](ch07.md)

-   **第8章 · 贝叶斯推断与 Gibbs 抽样**

    ---

    MLE 与 MAP、先验/似然/后验、共轭分布、后验预测分布、MCMC 算后验、Gibbs 抽样与二维正态实例。[S&M 6.1–6.2 ∪ D&S 6.2–6.4]

    [:octicons-arrow-right-24: 开始](ch08.md)

-   **第9章 · 蒙特卡罗积分与方差缩减**

    ---

    误差来源与确定性方法的高维失效、hit-or-miss 与样本均值估计器、控制变量与重要性抽样、MCMC 求积与舍入误差。[S&M 6.3–6.4 ∪ D&S 5]

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
