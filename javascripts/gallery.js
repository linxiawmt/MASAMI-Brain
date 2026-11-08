// 首页标记：给 body 加 is-home（全屏壁纸 + 玻璃拟态仅首页生效）
(function () {
  function setup() {

    // 首页视频背景：兜底自动播放（部分浏览器策略会拦截 muted autoplay）
    const bgVideo = document.querySelector(".home-bg-video");
    if (bgVideo) {
      bgVideo.muted = true;
      const tryPlay = () => {
        const p = bgVideo.play();
        if (p && p.catch) p.catch(() => {});
      };
      tryPlay();
      document.addEventListener("visibilitychange", () => {
        if (!document.hidden) tryPlay();
      });
      // 首次用户交互后重试
      ["click", "touchstart", "keydown"].forEach(function (ev) {
        document.addEventListener(ev, tryPlay, { once: true });
      });
    }
    const isHome = !!document.querySelector(".home-hero");
    if (isHome) {
      document.body.classList.add("is-home");
    } else {
      document.body.classList.remove("is-home");
    }
  }

  if (typeof document$ !== "undefined" && document$.subscribe) {
    document$.subscribe(setup);       // Material 即时导航，每次换页都校正
  } else {
    document.addEventListener("DOMContentLoaded", setup);
  }
})();
