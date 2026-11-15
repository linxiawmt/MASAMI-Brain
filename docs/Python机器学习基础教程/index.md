---
icon: material/book-education
---

# Python机器学习基础教程

Andreas C. Müller & Sarah Guido《Introduction to Machine Learning with Python》全书 **8 章**精读笔记：这是 scikit-learn 官方贡献者写的入门书，讲"怎么用 scikit-learn 把机器学习做对"——监督学习（kNN / 线性模型 / 朴素贝叶斯 / 决策树 / 随机森林 / 梯度提升 / SVM / 神经网络）→ 无监督与预处理（缩放 / PCA / NMF / 流形 / 聚类）→ 特征工程 → 交叉验证与网格搜索 → 评估指标 → Pipeline → 文本数据（词袋 / tf-idf / 主题建模）。

!!! note "运行环境"
    原书基于 Python 3.6 + scikit-learn 0.19，示例大量使用作者自己的绘图库 `mglearn`（`pip install mglearn`）。笔记中的代码按原书忠实呈现，**scikit-learn 1.x 下的 API 变化在正文用"版本差异"标注**（如 `load_boston` 在 1.2 起移除）。插图均截自原书。

!!! tip "怎么读这本书"
    和《机器学习实战》那种"从零手写算法"的路线不同，本书教的是**调用现成工具并把流程做规范**：训练/测试集划分 → fit/predict → 评估 → 调参 → 流水线。两本书正好互补，术语写法两边保持一致。

<div class="grid cards" markdown>

-   **第1章 · 引言**

    ---

    为何选机器学习、scikit-learn 生态、必要的库，以及第一个完整应用：鸢尾花三分类。

    [:octicons-arrow-right-24: 开始](ch01.md)

-   **第2章 · 监督学习**

    ---

    全书最长的一章：kNN、线性模型（回归/分类/正则化）、朴素贝叶斯、决策树与随机森林、梯度提升、SVM、神经网络，以及不确定度估计。

    [:octicons-arrow-right-24: 开始](ch02.md)

-   **第3章 · 无监督学习与预处理**

    ---

    缩放四种姿势、PCA 与白化、NMF 人脸分解、流形学习（t-SNE）、K-均值 / 凝聚聚类 / DBSCAN。

    [:octicons-arrow-right-24: 开始](ch03.md)

-   **第4章 · 数据表示与特征工程**

    ---

    分类变量 OneHotEncoder、分箱与交互特征、多项式特征、log 变换、三种特征选择、专家知识造特征（时间序列）。

    [:octicons-arrow-right-24: 开始](ch04.md)

-   **第5章 · 模型评估与改进**

    ---

    交叉验证五种划分、网格搜索嵌套 CV、两类错误与不平衡指标（ROC/AUC、PR 曲线）、多分类指标、回归指标。

    [:octicons-arrow-right-24: 开始](ch05.md)

-   **第6章 · 算法链与管道**

    ---

    Pipeline：把缩放器和模型绑成一个估计器，杜绝"缩放泄漏"；管道 + 网格搜索一次搜全流程超参。

    [:octicons-arrow-right-24: 开始](ch06.md)

-   **第7章 · 处理文本数据**

    ---

    词袋、停用词、tf-idf、n 元分词、模型系数解读、高级分词与词干提取、LDA 主题建模。

    [:octicons-arrow-right-24: 开始](ch07.md)

-   **第8章 · 全书总结**

    ---

    机器学习方法论：面对新问题怎么下手、原型到生产、测试系统、自己写估计器、进阶资源。

    [:octicons-arrow-right-24: 开始](ch08.md)

</div>
