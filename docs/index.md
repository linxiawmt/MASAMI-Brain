---
hide:
  - toc
  - navigation
---

<div class="home-bg-fixed" style="background-image:url('images/bg/ellen.jpg')">
  <video class="home-bg-video" autoplay muted loop playsinline preload="metadata" poster="images/bg/miku-poster.jpg">
    <source src="images/bg/miku-loop.mp4" type="video/mp4">
  </video>
</div>
<div class="home-bg-glow"></div>

<div class="home-wrap">

  <!-- HERO：整行对称，左侧自述 + 右侧"手头在啃" -->
  <section class="glass bento hero reveal" style="--i:0">
    <div class="hero-left">
      <div class="home-kicker">记在网上的学习笔记本</div>
      <div class="home-title">MASAMI 的外置大脑</div>
      <div class="home-sub">啃过的书、推导过的公式、跑通过的代码都堆在这里——先给自己回头翻，顺手让同路的人少踩几个坑。</div>
      <div class="home-btns">
        <a class="gbtn" href="ml/">去翻最新一本</a>
        <a class="gbtn ghost" href="https://github.com/linxiawmt/MASAMI-Brain" target="_blank" rel="noopener">GitHub</a>
      </div>
    </div>
    <aside class="hero-now">
      <div class="now-kicker">NOW · 手头在啃</div>
      <div class="now-item">
        <div class="now-name">《机器学习实战》</div>
        <div class="now-note">15 章过了 9 章，剩下的在路上</div>
        <div class="now-bar"><i style="--w:60%"></i></div>
      </div>
      <div class="now-item">
        <div class="now-name">《传热学》</div>
        <div class="now-note">14 章完结留档，过阵子回头复习</div>
        <div class="now-bar done"><i style="--w:100%"></i></div>
      </div>
    </aside>
  </section>

  <!-- 数字瓷片：四个等宽，数字滚动 -->
  <div class="bento stat reveal" style="--i:1">
    <div class="stat-num" data-count="8">0</div>
    <div class="stat-label">书架板块</div>
    <div class="stat-sub">SECTIONS</div>
  </div>
  <div class="bento stat reveal" style="--i:2">
    <div class="stat-num" data-count="93">0</div>
    <div class="stat-label">章笔记在线</div>
    <div class="stat-sub">CHAPTERS</div>
  </div>
  <div class="bento stat reveal" style="--i:3">
    <div class="stat-num" data-count="3">0</div>
    <div class="stat-label">套代码栈</div>
    <div class="stat-sub">PYTHON · MATLAB · PYTORCH</div>
  </div>
  <div class="bento stat reveal" style="--i:4">
    <div class="stat-num" data-count="2">0</div>
    <div class="stat-label">副深浅主题</div>
    <div class="stat-sub">LIGHT · DARK</div>
  </div>

</div>

  <!-- 书架：九张等宽瓷片，3×3 对称（markdown 容器必须在顶层） -->
  <div class="sec-head reveal" style="--i:5">
    <span class="sec-kicker">PROJECTS</span>
    <span class="sec-title">书架上现在有九样东西</span>
  </div>

<div class="home-grid-cards" markdown>
<div class="grid cards" markdown>

-   :fontawesome-solid-eye:{ .lg .middle } **深度学习视觉系统**

    ---

    Elgendy 那本 10 章啃完，Keras 示例全翻成了 PyTorch——能跑的我都自己跑过。

    [:octicons-arrow-right-24: 去看](dl4vs/index.md)

-   :fontawesome-solid-calculator:{ .lg .middle } **数值分析方法**

    ---

    13 章理论为主，推导不跳步，MATLAB 例程直接内联在页面上。期末复习靠它续命。

    [:octicons-arrow-right-24: 去看](numana/index.md)

-   :fontawesome-solid-wave-square:{ .lg .middle } **计算傅里叶光学**

    ---

    Voelz (SPIE TT89) 9 章，衍射、透镜、成像仿真，MATLAB 和 Python 双份代码对照写。

    [:octicons-arrow-right-24: 去看](holo/fourier/index.md)

-   :fontawesome-solid-wind:{ .lg .middle } **光波传播数值模拟**

    ---

    Schmidt (SPIE PM199)，采样判据、相位屏这些硬骨头都在；书缺的第 3 章留了占位。

    [:octicons-arrow-right-24: 去看](holo/waveprop/index.md)

-   :fontawesome-solid-camera-retro:{ .lg .middle } **光学扫描全息**

    ---

    Poon (Springer) 5 章，OSH 从原理到进展，公式全推了一遍。

    [:octicons-arrow-right-24: 去看](holo/scanholo/index.md)

-   :fontawesome-solid-images:{ .lg .middle } **OpenCV 轻松入门**

    ---

    李立宗那本 19 章，例子里能跑的全实跑过——页面截图是我机器上跑出来的，不是翻拍。

    [:octicons-arrow-right-24: 去看](opencv/index.md)

-   :fontawesome-solid-fire:{ .lg .middle } **传热学**

    ---

    Incropera 圣经 14 章，例题按"已知 → 求 → 分析 → 点评"完整拆解。

    [:octicons-arrow-right-24: 去看](heat/index.md)

-   :fontawesome-solid-robot:{ .lg .middle } **机器学习实战**

    ---

    Harrington 的 MLIA，9/15 章上线，代码全从 Python 2 改成了能直接跑的 Python 3。

    [:octicons-arrow-right-24: 去看](ml/index.md)

-   :fontawesome-solid-compass:{ .lg .middle } **学习资源导航**

    ---

    顺手收藏的 GitHub 好东西：浙大专属、CS 自学、AI……和课程笔记分开放。

    [:octicons-arrow-right-24: 去看](resources/index.md)

</div>
</div>

<div class="home-tail">

  <!-- 底部两卡：等宽对称 -->
  <section class="glass bento bcard reveal" style="--i:15">
    <div class="hc-title">怎么用</div>
    <div class="hc-item">代码书：页面笔记 + 能跑的代码 + 实跑截图，三者永远对得上</div>
    <div class="hc-item">理论书：公式全推导，原书跳步的地方我补了，笔误直接在页面里标</div>
    <div class="hc-item">顶部标签换栏目 · 左侧翻章节 · 右上角放大镜搜全文 · 深浅色随手切</div>
  </section>

  <section class="glass bento bcard reveal" style="--i:16">
    <div class="hc-title">关于这个站</div>
    <div class="hc-item">MkDocs Material 驱动，界面 CSS 自己手搓，背景是 Wallpaper Engine 里那张 MIKU</div>
    <div class="hc-item">笔记忠于原书，发现原书笔误就直接标出来改掉，不藏着</div>
    <div class="hc-item">更新不定期，但一直在写——最新进度看各栏目入口页</div>
  </section>

</div>
