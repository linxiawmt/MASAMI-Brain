// 首页标记：给 body 加 is-home（全屏壁纸 + 玻璃拟态仅首页生效）
(function () {
  function setup() {
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
