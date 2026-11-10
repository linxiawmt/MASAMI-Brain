# 附录 A/B 生成均匀随机数 · Perron–Frobenius 定理

> **对应教材**：Shonkwiler & Mendivil《Explorations in Monte Carlo Methods》附录 A（pp.245–254）、附录 B（pp.255–258）。附录 A 是 §1.3 随机数生成的深水区（多存储值发生器、有限域构造、Mersenne Twister、两种统计检验）；附录 B 是 §3.1 不变分布收敛性（式 (3.6)）背后定理的完整证明。附录 C/D 的结论已并入 ch04（§4.6.4 拓展框、§4.3 Donsker 提及），附录 E（项目清单）从略。
> **补充**：Dunn & Shultis 的 RNG 章节与其检验方法可与 A.3 对照阅读。

## 附录 A 生成均匀随机数（Generating Uniform Random Numbers）

### A.1 多存储值随机数生成（Multiple Stored Value RNG）

§1.3.3 的 LCRNG 周期至多 \(m\)——**单种子**的宿命：种子一重复全序列重复。用**两个以上存储值**可让周期远超 \(m\)。

**Fibonacci 型发生器**：双种子 + 参数 \(a_1,a_2,c\) 与模 \(m\)，\(x_1=seed_1,\ x_0=seed_2\)：

\[
x_{n+1}=(a_1x_n+a_2x_{n-1}+c)\bmod m\qquad(\text{常取 }c=0).
\]

名字致敬 Fibonacci 数列 \(x_n=x_{n-1}+x_{n-2}\)。推广成**存储值表**：好成绩的发生器

\[
x_n=(x_{n-24}+x_{n-55})\bmod m,\qquad n\ge55,
\]

（两个乘子均为 1；需预填 \(x_0,\dots,x_{54}\)。）书 Fig A.1 给出 Marsaglia 的实现（C 代码），关键一行是"next = 某旧值 − 另一旧值"的 Fibonacci 形态，新值回填进表。

**周期为何能爆表**：状态不再是单个 \(x_n\) 而是**整张表**（如 55 个值）——状态空间 \(m^{55}\)，周期可天文。

### A.1.2 有限域 RNG（Finite Field RNG）

用有限域理论构造**双存储值、周期恰 \(m^2-1\)** 的发生器（书以 \(m=5\) 小例演示）。构造地在域

\[
G_{25}=\{a_0+a_1x:\ a_0,a_1\in\mathbb Z_5,\ x^2=3\}
\]

（\(\mathbb Z_5\) 上模 \(x^2+2\) 的多项式；\(x^2+2\) 在 \(\mathbb Z_5\) 上不可约故确为域）。要找**本原元**（primitive element）\(\beta=a_0+a_1x\)，使 \(\beta\) 的幂跑遍 \(G_{25}\) 的全部非零元（共 \(25-1=24\) 个）。书试 \(\beta=2+3x\)：验证它是 \(z^2-4z-3\) 的根（\(\beta^2-4\beta-3=(4+12x+9x^2)-(8+12x)-3=3+4x^2=3+4\times3=0\)），于是 \(\beta^2=4\beta+3\)；列出 \(\beta^1\dots\beta^{24}\) 的表确认 24 个幂遍历非零元。

把"系数求和"映射回 \(\mathbb Z_5\)（如 \(4+2\beta\to4+2=1\bmod5\)）即得输出序列。更优雅的等价物：由 \(\beta^2=4\beta+3\) 直接写线性递推（对连续三个幂 \(\beta^{k+2}=4\beta^{k+1}+3\beta^k\) 系数求和）：

\[
y_n=(4y_{n-1}+3y_{n-2})\bmod5,
\]

初值 \(y_0=4,y_1=1\) 给出 24 长序列 `1 2 1 0 3 2 2 4 2 0 1 4 4 3 4 0 2 3 3 1 3 0 4 1`。周期 \(m^2-1=24\)（不可达 25，脚注习题：素模 \(m\) 下 \(x\mapsto ax+c\) 必有不动点）。

### A.2 Mersenne Twister

Maple/MATLAB/Python 等的标配伪随机数发生器：\(w\) 位字长下生成 \([0,2^w-1]\) 均匀整数；**快、过大多数统计检验、无专利、周期 = 梅森素数 \(2^{19937}-1\)**（6002 位十进制）。

内部计算几乎全是**位运算**：由种子（\(x_0\)=seed）递推 624 个 \(w\) 位整数的状态表：

\[
x_i=f\ast(x_{i-1}\oplus(x_{i-1}\gg(w-2)))+i,
\]

（\(\oplus\) 异或、\(\gg\) 右移；\(f=1812433253\)。）输出时每个表值先过"tempering（回火）"函数；623 个用完后，最后一个经"twist（扭转）"生成新表。Python 用户：就是发行目录里的 `random.py`。

### A.3 非随机性检验（Testing for Non-randomness）

> 检验只能证明序列**非**随机，永远不能证明它随机。更多检验见 [Knu81]。

**χ² 检验**：结果分成 \(k\) 类，假设概率 \(p_i\)（\(\sum p_i=1\)）；\(N\) 次试验中类 \(i\) 的期望 \(p_iN\)、实际 \(Y_i\)：

\[
V=\sum_{i=1}^{k}\frac{(Y_i-p_iN)^2}{p_iN}
\]

是自由度 \(k-1\) 的 χ² 随机变量。书例：20 次掷双硬币 HH=2, HT=5, TH=3, TT=10 → \(V=7\)；χ² 表（\(\nu=3\)）给 \(\Pr(V\le7)\approx0.906\)，即 \(P(V\ge7)=9.4\%\)——90% 水平下不过关、95% 下过关。

**Kolmogorov–Smirnov 检验**：检验样本 \(x_1,\dots,x_N\) 是否来自 cdf \(F\)。算**经验 cdf** \(\hat F_N(t)=\#\{i:x_i\le t\}/N\)；Glivenko–Cantelli 定理保证真来自 \(F\) 时 \(\sup_t|\hat F_N-F|\to0\)。统计量

\[
D_N=\sup_t\big|\hat F_N(t)-F(t)\big|
\]

在独立性假设下分布已知、可查表（书 Table A.1：\(n=12\)、\(\alpha=0.05\) 的门限 0.38）。书例：12 个数据点的 \(D_{12}=0.239\) → 过关；LCRNG \(x_{n+1}=1+3x_n\bmod p\) 三个模数的对照：\(p=11\)（周期 5，1000 样本绕 200 圈）\(D=0.345\) 远超门限 \(1.36/\sqrt{1000}=0.043\)——挂；\(p=229\)：\(D=0.146\)——挂；\(p=104729\)：\(D=0.024\)——过（但这只说明"不能断定它不随机"）。

## 附录 B Perron–Frobenius 定理

> 本附录证明的版本特化为**行随机矩阵**（直接适用于马尔可夫链转移矩阵）；一般版本只要求非负。证明路线：振荡递缩 → 列收敛 → 特征值分离 → Jordan 形式收尾。

!!! note "定理 B.1（Perron–Frobenius，行随机版）"
    设 \(P\) 是 \(n\times n\) 行随机矩阵且 \(P_{i,j}\ge\delta>0\) 对一切 \(i,j\)。则：

    1. \(P\) 有特征值 1，对应**规范化左特征向量 \(\pi^T\) 严格为正**；
    2. 其余特征值的模严格小于 1；
    3. \(\displaystyle\lim_{n\to\infty}P^n=\mathbb 1\pi^T\)（\(\mathbb 1\)=全 1 列向量），从而对一切概率向量 \(\alpha\)：\(\lim\alpha P^n=\pi\)。

    \(\pi\) 唯一——这就是 \(\pi P=\pi\) 存在唯一且 \(p_t=p_0P^t\to\pi\)（式 (3.5)(3.6)）的依据。

**证明骨架**（书 B.1）：

- **不变子空间**：\(P\mathbb 1=\mathbb 1\)；\(\langle f,\mathbb 1\rangle=0\Rightarrow\langle P^Tf,\mathbb 1\rangle=\langle f,P\mathbb 1\rangle=0\)，故 \(\mathcal S=\{v:\langle v,\mathbb 1\rangle=0\}\) 是 \(P^T\) 的不变子空间。
- **振荡压缩（Prop B.2）**：定义 \(\operatorname{osc}(f)=\max f-\min f\)。条件 \(P_{ij}\ge\delta\) 下，把 \((Pf)_i\) 拆成"\(\delta\) 份最小值 + 其余"：

  \[
  (Pf)_i\le\delta f_a+(1-\delta)f_b\quad\text{且}\quad (Pf)_i\ge(1-\delta)f_a+\delta f_b
  \quad\Longrightarrow\quad
  \operatorname{osc}(Pf)\le(1-2\delta)\operatorname{osc}(f).
  \]

  （\(\delta\) 的保证：不可约有限链的某个 \(P\) 幂满足它——条件并不苛刻。）
- **列收敛（Prop B.3）**：记 \(c_j^n\) 为 \(P^n\) 的第 \(j\) 列，则 \(\operatorname{osc}(c_j^{n+1})\le(1-2\delta)\operatorname{osc}(c_j^n)\to0\)，故 \(c_j^n\to\pi_j\mathbb 1\)；几何速率 \(C(1-2\delta)^n\)。\(P^n\to\mathbb 1\pi^T\)；由 \((\mathbb 1\pi^T)P=P^{n+1}\) 的极限得 \(\pi^TP=\pi^T\)。
- **特征值分离（Prop B.4）**：对 \(f\in\mathcal S\)，\((\hat{P^T})^nf\to\pi\langle\mathbb 1^Tf\rangle=0\)——\(\hat{P^T}\)（\(P^T\) 在 \(\mathcal S\) 上的限制）一切特征值模 <1；\(P^T\) 的特征值集合 = {1} ∪ 这些，故 \(P\) 只有 1 的模为 1。
- **推论 B.5**：\(P\) 不可约且非周期（不必正）时，某幂 \(P^n\) 为正矩阵，对 \(P^n\) 用定理：唯一正随机行向量 \(\pi\) 使 \(\pi P=\pi\)（用 \(P^3\) 的情形做唯一性桥接：\(\phi_2P=\phi_1P^2=\pi P^3=\pi\) 与 \(\phi_2P^3=\pi\) 联立逼出 \(\pi P=\pi\)）。
- **观察（Prop B.6）**：Jordan 形式 \(P=C^{-1}JC\) 下，\(J^k\to\) 单 1 矩阵 \(L\)，且 \(C^{-1}\) 的第一列恰是 \(\mathbb 1\)（由 \(C C^{-1}=I\) 的点积约束 + 唯一解），于是 \(\lim P^k=C^{-1}LC=\mathbb 1\pi\)——外积形式收尾。

!!! note "与第 3/4/5 章的接缝"
    - §3.1.2 "正则链收敛到不变分布" = 定理 B.1 + 推论 B.5（不可约 + 非周期 ⇒ 某幂为正）。
    - §3.2.2 Metropolis 链正则性的依据：自环 ⇒ 非周期、提案不可约 ⇒ 链不可约。
    - §4.7.2 逃逸速率 \(\omega\hat P=\lambda\omega\)、\(\lambda<1\)（式 (4.33)）= Prop B.4 的瞬态版（行和 <1 的子阵）。
    - ch01 习题 2 的谱隙（次特征值控制收敛速率）= 证明中 \(1-2\delta\) 几何速率的精细版。

## 附录 C/D/E 去向说明

- **附录 C（相关投资的 Kelly 分配）**：结论与四项增长率公式已并入 [ch04 §4.6.4 拓展框](ch04.md)；书 C.2 的 GA 代码属工程实现，从略。
- **附录 D（Donsker 定理）**：随机游走的连续插值 \(W(n,t)\)（式 (D.1)）弱收敛到 Wiener 过程——[ch04 §4.3](ch04.md) 已引（"布朗运动常返性由 Pólya + Donsker 导出"）；函数式收敛陈述见原书。
- **附录 E（Projects）**：课程项目选题清单（TSP、clique、蛋白质结构、Monopoly、put-call parity 验证、人工生态、Mini-hearts、分形逆问题等），不设笔记页；做课程设计时直接翻原书 pp.267–270。
