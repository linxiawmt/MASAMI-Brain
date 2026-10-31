// 首页画廊灯箱：点击卡片看大图，点击任意处 / Esc 关闭
(function () {
  function setup() {
    if (document.querySelector(".lb")) return;
    const imgs = document.querySelectorAll(".g-strip .g-card img");
    if (!imgs.length) return;

    const lb = document.createElement("div");
    lb.className = "lb";
    lb.innerHTML = '<img alt="preview">';
    document.body.appendChild(lb);
    const big = lb.querySelector("img");

    function open(src) {
      big.src = src;
      lb.classList.add("open");
    }
    function close() {
      lb.classList.remove("open");
    }

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

  if (typeof document$ !== "undefined" && document$.subscribe) {
    document$.subscribe(setup);       // Material 即时导航
  } else {
    document.addEventListener("DOMContentLoaded", setup);
  }
})();
