---
icon: material/robot-industrial
---

# 机器学习实战

Peter Harrington《Machine Learning in Action》全书 **15 章**精读笔记：分类（kNN / 决策树 / 朴素贝叶斯 / Logistic / SVM / AdaBoost）→ 回归（线性 / 树回归）→ 无监督（K-均值 / Apriori / FP-growth）→ 工具（PCA / SVD / MapReduce）。

!!! note "运行环境"
    conda 环境 `machinelearn`（Python 3.10 + NumPy + Matplotlib）。原书代码是 Python 2 风格，笔记中已全部改写为可直接运行的 Python 3，差异点在正文标注。程序运行结果插图均来自本机实跑（标注"代码实跑"），概念示意图截自原书。

<div class="grid cards" markdown>

-   **第1章 · 机器学习基础**

    ---

    机器学习是什么、关键术语、任务类型与选型，NumPy 函数库基础。

    [:octicons-arrow-right-24: 开始](ch01.md)

-   **第2章 · k-近邻算法**

    ---

    距离投票分类：约会网站配对、手写识别，归一化与误差率评估。

    [:octicons-arrow-right-24: 开始](ch02.md)

-   **第3章 · 决策树**

    ---

    信息增益划分、递归建树、Matplotlib 绘树、隐形眼镜实例与 pickle 存储。

    [:octicons-arrow-right-24: 开始](ch03.md)

-   **第4章 · 朴素贝叶斯**

    ---

    条件概率分类：词向量、伯努利与词袋模型、垃圾邮件过滤与 RSS 地域倾向。

    [:octicons-arrow-right-24: 开始](ch04.md)

-   **第5章 · Logistic回归**

    ---

    Sigmoid 分类、梯度上升与随机梯度上升找最佳系数、疝气病马死亡率预测。

    [:octicons-arrow-right-24: 开始](ch05.md)

-   **第6章 · 支持向量机**

    ---

    最大间隔、SMO 优化、完整 Platt SMO、径向基核函数与手写识别回顾。

    [:octicons-arrow-right-24: 开始](ch06.md)

-   **第7章 · 利用AdaBoost元算法提高分类性能**

    ---

    bagging 与 boosting、单层决策树弱分类器、ROC 曲线与非均衡分类代价。

    [:octicons-arrow-right-24: 开始](ch07.md)

-   **第8章 · 预测数值型数据：回归**

    ---

    线性回归、局部加权回归（LWLR）、岭回归与 lasso、鲍鱼年龄与乐高价格。

    [:octicons-arrow-right-24: 开始](ch08.md)

-   **第9章 · 树回归**

    ---

    CART 建树、预剪枝与后剪枝、模型树、树回归 vs 标准回归与 Tkinter GUI。

    [:octicons-arrow-right-24: 开始](ch09.md)

-   **第10章 · 利用K-均值聚类算法对未标注数据分组**

    ---

    K-均值流程、SSE 与后处理、二分 K-均值、地图坐标聚类。

    [:octicons-arrow-right-24: 开始](ch10.md)

-   **第11章 · 使用Apriori算法进行关联分析**

    ---

    频繁项集与支持度、Apriori 原理剪枝、关联规则挖掘、国会投票与毒蘑菇。

    [:octicons-arrow-right-24: 开始](ch11.md)

-   **第12章 · 使用FP-growth算法来高效发现频繁项集**

    ---

    FP 树编码、条件模式基与条件树、Twitter 共现词与点击流挖掘。

    [:octicons-arrow-right-24: 开始](ch12.md)

-   **第13章 · 利用PCA来简化数据**

    ---

    降维技术、协方差矩阵特征分解、方差百分比取舍、半导体制造数据降维。

    [:octicons-arrow-right-24: 开始](ch13.md)

-   **第14章 · 利用SVD简化数据**

    ---

    隐性语义索引、协同过滤与相似度、SVD 增强推荐、图像压缩。

    [:octicons-arrow-right-24: 开始](ch14.md)

-   **第15章 · 大数据与MapReduce**

    ---

    Hadoop 流、mrjob 自动化 MapReduce、分布式 SVM 的 Pegasos 算法。

    [:octicons-arrow-right-24: 开始](ch15.md)

</div>
