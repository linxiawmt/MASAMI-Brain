// 首页动效合集：is-home 标记 / 视频壁纸兜底播放 / 入场 reveal /
// 数字滚动 count-up / 卡片鼠标聚光。兼容 Material instant 导航。
(function () {
  "use strict";

  var reducedMotion = window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function easeOutCubic(t) { return 1 - Math.pow(1 - t, 3); }

  // 数字滚动：0 → data-count，仅首次
  function countUp(el) {
    if (el.dataset.counted) return;
    el.dataset.counted = "1";
    var target = parseInt(el.dataset.count, 10) || 0;
    if (reducedMotion || !window.requestAnimationFrame) {
      el.textContent = String(target);
      return;
    }
    var dur = 950, t0 = null;
    function step(ts) {
      if (t0 === null) t0 = ts;
      var p = Math.min((ts - t0) / dur, 1);
      el.textContent = String(Math.round(easeOutCubic(p) * target));
      if (p < 1) window.requestAnimationFrame(step);
    }
    window.requestAnimationFrame(step);
  }

  function setup() {
    var isHome = !!document.querySelector(".home-hero, .bento.hero");
    document.body.classList.toggle("is-home", isHome);

    // JS 就绪标记：入场动画样式只在此时生效（无 JS 也能正常阅读）
    document.body.classList.add("fx-ready");

    // 视频壁纸：兜底自动播放
    var bgVideo = document.querySelector(".home-bg-video");
    if (bgVideo) {
      bgVideo.muted = true;
      var tryPlay = function () {
        var p = bgVideo.play();
        if (p && p.catch) p.catch(function () {});
      };
      tryPlay();
      document.addEventListener("visibilitychange", function () {
        if (!document.hidden) tryPlay();
      });
      ["click", "touchstart", "keydown"].forEach(function (ev) {
        document.addEventListener(ev, tryPlay, { once: true });
      });
    }

    if (!isHome) return;

    // 书架瓷片统一编入场序号（每张错开 60ms）
    var cards = document.querySelectorAll(
      ".home-grid-cards .grid.cards > ul > li, .home-grid-cards .grid.cards > li"
    );
    cards.forEach(function (li, i) { li.style.setProperty("--i", i); });

    // 入场 reveal：进视口加 .in；数字块同时触发 count-up
    var targets = [];
    document.querySelectorAll(".reveal").forEach(function (el) { targets.push(el); });
    cards.forEach(function (el) { targets.push(el); });

    if (reducedMotion || typeof IntersectionObserver === "undefined") {
      targets.forEach(function (el) {
        el.classList.add("in");
        el.querySelectorAll(".stat-num[data-count]").forEach(countUp);
      });
      return;
    }

    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        if (!entry.isIntersecting) return;
        entry.target.classList.add("in");
        entry.target.querySelectorAll(".stat-num[data-count]").forEach(countUp);
        io.unobserve(entry.target);
      });
    }, { threshold: 0.12, rootMargin: "0px 0px -8% 0px" });
    targets.forEach(function (el) { io.observe(el); });

    // 卡片鼠标聚光：把光标坐标写进 CSS 变量
    document.querySelectorAll(".bento").forEach(function (card) {
      card.addEventListener("mousemove", function (e) {
        var r = card.getBoundingClientRect();
        card.style.setProperty("--mx", (e.clientX - r.left) + "px");
        card.style.setProperty("--my", (e.clientY - r.top) + "px");
      });
    });
  }

  if (typeof document$ !== "undefined" && document$.subscribe) {
    document$.subscribe(setup);   // Material 即时导航，每次换页都校正
  } else {
    document.addEventListener("DOMContentLoaded", setup);
  }
})();
