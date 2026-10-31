// 首页画廊灯箱 + Wallhaven 唯美壁纸轮换 + 首页全页壁纸底
(function () {
  function setup() {
    if (document.querySelector(".lb")) return;

    // ---- 灯箱 ----
    const imgs = document.querySelectorAll(".g-strip .g-card img");
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

    // ---- 主页：唯美壁纸轮换（Wallhaven 高分精选）----
    const hero = document.querySelector(".hero");
    if (hero) {
      document.body.classList.add("is-home");
      const POOL = [
        "images/wallhaven/a1.jpg",
        "images/wallhaven/a2.jpg",
        "images/wallhaven/a3.jpg",
        "images/wallhaven/a4.jpg",
        "images/wallhaven/a5.jpg",
        "images/wallhaven/a6.jpg",
        "images/wallhaven/a7.jpg",
      ];
      const layers = [...hero.querySelectorAll(".hero-bg")];
      const dotsBox = hero.querySelector(".hero-dots");
      const dots = POOL.map(function (_, i) {
        const d = document.createElement("span");
        d.className = "dot" + (i === 0 ? " on" : "");
        dotsBox.appendChild(d);
        return d;
      });
      let cur = 0;
      const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
      function show(i) {
        if (i === cur) return;
        const nxt = layers.find(function (l) { return !l.classList.contains("on"); });
        const pre = new Image();
        pre.src = POOL[i];                       // 先预加载，避免切换闪空
        nxt.style.backgroundImage = "url('" + POOL[i] + "')";
        layers.forEach(function (l) { l.classList.toggle("on", l === nxt); });
        dots.forEach(function (d, j) { d.classList.toggle("on", j === i); });
        cur = i;
      }
      if (!reduce && POOL.length > 1) {
        setInterval(function () { show((cur + 1) % POOL.length); }, 9000);
      }
    }
  }

  if (typeof document$ !== "undefined" && document$.subscribe) {
    document$.subscribe(setup);       // Material 即时导航
  } else {
    document.addEventListener("DOMContentLoaded", setup);
  }
})();
