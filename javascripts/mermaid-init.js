// Mermaid 本地渲染（v10）：fence_div_format 输出 <div class="mermaid">源码</div>。
// 本脚本是唯一渲染者：显式关闭 startOnLoad，随即时导航重绘，跟随深浅主题。
// 渲染失败时把源码放回并插入可见提示，避免黑盒吞错。
(function () {
  var lastScheme = null;
  function renderAll(force) {
    if (!window.mermaid) return;
    var schemeEl = document.body || document.documentElement;
    var dark = schemeEl.getAttribute("data-md-color-scheme") === "slate";
    if (dark === lastScheme && !force) return;
    lastScheme = dark;
    mermaid.initialize({
      startOnLoad: false,
      securityLevel: "loose",
      theme: dark ? "dark" : "neutral",
      flowchart: { curve: "basis", htmlLabels: true }
    });
    document.querySelectorAll("div.mermaid").forEach(function (el) {
      var src = el.dataset.mmdSrc || (el.dataset.mmdSrc = el.textContent);
      if (force) {
        el.dataset.mmdDone = "";
        el.textContent = src;
      }
      if (el.dataset.mmdDone) return;
      el.dataset.mmdDone = "1";
      if (!src.trim()) return;
      var id = "mmd-" + Math.random().toString(36).slice(2);
      mermaid.render(id, src).then(function (res) {
        el.innerHTML = res.svg;
      }).catch(function (e) {
        var msg = e && (e.message || e.str) || String(e);
        el.textContent = src;
        var warn = document.createElement("div");
        warn.style.cssText = "color:#e07b54;font-size:12px;padding:4px";
        warn.textContent = "[mermaid] " + msg.slice(0, 200);
        el.appendChild(warn);
      });
    });
  }
  if (window.mermaid) {
    mermaid.initialize({ startOnLoad: false });
  }
  document$.subscribe(function () { renderAll(true); });
  // 深浅色切换时让流程图跟随主题重渲染（scheme 属性挂在 body 上）
  window.addEventListener("DOMContentLoaded", function () {
    new MutationObserver(function () { renderAll(false); }).observe(
      document.body, { attributes: true, attributeFilter: ["data-md-color-scheme"] }
    );
  });
})();
