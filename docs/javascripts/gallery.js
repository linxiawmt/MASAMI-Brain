// 首页画廊灯箱：点击卡片看大图，点击任意处 / Esc 关闭
(function () {
  function setup() {
    if (document.querySelector(".lb")) return;
    const imgs = document.querySelectorAll(".g-strip .g-card img");

    // 灯箱
    if (imgs.length) {
      const lb = document.createElement("div");
      lb.className = "lb";
      lb.innerHTML = '<img alt="preview">';
      document.body.appendChild(lb);
      const big = lb.querySelector("img");
      const open = (src) => { big.src = src; lb.classList.add("open"); };
      const close = () => lb.classList.remove("open");
      imgs.forEach(function (el) {
        el.addEventListener("click", function () {
          open(el.currentSrc || el.src);
        });
      });
      lb.addEventListener("click", close);
      document.addEventListener("keydown", function (e) {
        if (e.key === "Escape") close();
      });
    }

    // 主页：Wallhaven 壁纸轮换 + 全页壁纸底
    const heroBgs = document.querySelectorAll(".hero .hero-bg");
    if (heroBgs.length) {
      document.body.classList.add("is-home");
      let cur = 0;
      heroBgs[0].classList.add("on");
      if (heroBgs.length > 1 && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
        setInterval(function () {
          heroBgs[cur].classList.remove("on");
          cur = (cur + 1) % heroBgs.length;
          heroBgs[cur].classList.add("on");
        }, 9000);
      }
    }
  }

  if (typeof document$ !== "undefined" && document$.subscribe) {
    document$.subscribe(setup);       // Material 即时导航
  } else {
    document.addEventListener("DOMContentLoaded", setup);
  }
})();
