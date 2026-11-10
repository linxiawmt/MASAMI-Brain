# 随机模拟方法及应用 · 笔记总览

> **骨架教材**：Shonkwiler & Mendivil《Explorations in Monte Carlo Methods》(2nd ed, Springer UTM, 2024)。
> **补充教材**：Dunn & Shultis《Exploring Monte Carlo Methods》——只在对应位置以"补充"框出现。
> **本节承诺线性**：笔记第 N 章 = 教材第 N 章，小节编号、顺序跟书走，可以书、笔记两头对着读。书里公式编号（如 式(1.25)）在笔记中原样引用。

## 章节对照（线性）

| 笔记 | 教材章 | 核心内容 | 配套代码 |
| --- | --- | --- | --- |
| [第 1 章](ch01.md) | Ch 1 导论（pp.1–46） | 蒲丰投针、直方图与样本路径、概率基础全推导（期望/方差/切比雪夫/贝叶斯/联合分布）、随机数生成（LCRNG）、hit-or-miss、优惠券收集、洗牌 | [ch01.ipynb](nb/ch01.ipynb) |
| [第 2 章](ch02.md) | Ch 2 概率分布及其用途（pp.47–102） | CDF 反演、Walker alias、二项/泊松/指数/柯西、中心极限定理与正态（Box–Muller 证明）、Gibrat 定律、拒绝采样与 beta、复合分布与 gamma、联合分布抽样 | [ch02.ipynb](nb/ch02.ipynb) |
| [第 3 章](ch03.md) | Ch 3 MCMC（pp.103–140） | 离散马尔可夫链、Metropolis 算法及正确性、遍历定理、统计力学与 Ising、Metropolis–Hastings、计数问题、coupling from the past | [ch03.ipynb](nb/ch03.ipynb) |
| [第 4 章](ch04.md) | Ch 4 随机游走（pp.141–184） | 1d 游走与常返性、扩散、布朗运动、期权定价、自回避行走、赌徒破产（差分方程全解）、Kelly 准则、电网与基本矩阵、kMC | [ch04.ipynb](nb/ch04.ipynb) |
| [第 5 章](ch05.md) | Ch 5 蒙特卡罗优化（pp.185–207） | 模拟退火（Hajek 冷却定理）、TSP 与路径部分反转、遗传算法、Shekel 函数、permanent 问题 | [ch05.ipynb](nb/ch05.ipynb) |
| [第 6 章](ch06.md) | Ch 6 再谈 MCMC（pp.209–244） | 贝叶斯推断（MAP/后验预测）、PyMC、Gibbs 抽样、MC 积分/求积、方差缩减（控制变量、重要性采样）、MCMC 求积、舍入误差 | [ch06.ipynb](nb/ch06.ipynb) |
| [附录 A/B](appendix.md) | App A+B（pp.245–258） | 生成均匀随机数（Fibonacci 型、有限域、Mersenne Twister、χ²/KS 检验）；Perron–Frobenius 定理证明 | — |

## 状态

| 章 | 笔记 | 代码 | 状态 |
| --- | --- | --- | --- |
| ch01 | 6.3 万字节 | ✅ 实跑 | ✅ 完成 |
| ch02 | 5.1 万字节 | ✅ 实跑 | ✅ 完成 |
| ch03 | — | — | ⬜ 未开工 |
| ch04 | — | — | ⬜ 未开工 |
| ch05 | — | — | ⬜ 未开工 |
| ch06 | — | — | ⬜ 未开工 |
| 附录 | — | — | ⬜ 未开工 |

!!! note "旧版去哪了"
    2026-11-10 之前本节曾有一版 11 章的主题式笔记（按肖柳青中文教材章序组织）。为配合英文教材线性学习，已整体归档到本地 `_mc_work/mc_old_archive/`，不再进站。
