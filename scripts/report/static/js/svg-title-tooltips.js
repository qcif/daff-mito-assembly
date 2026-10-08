// Promotes an inline SVG's native <title> tooltips to the page's
// Bootstrap tooltip styling, so graph nodes look like every other
// hover tooltip on the page instead of the plain browser popover.
// The standalone SVG file (task 47) keeps its <title> children —
// this only touches the copy embedded in report.html.
// `customClass` (optional) sets Bootstrap 5.0.2's `data-bs-custom-class`
// attribute, so a caller can style that tooltip's `.tooltip-inner`
// differently -- e.g. the genome map's multi-line tooltip text, which
// needs `white-space: pre-line` instead of the default single-line
// collapse (task 54 §5.5). The graph call omits it and is unaffected.
function promoteSvgTitleTooltips(containerId, itemSelector, customClass) {
  var container = document.getElementById(containerId);
  if (!container) return;
  var items = container.querySelectorAll(itemSelector);
  for (var i = 0; i < items.length; i++) {
    var el = items[i];
    var titleEl = null;
    for (var j = 0; j < el.children.length; j++) {
      if (el.children[j].tagName === 'title') {
        titleEl = el.children[j];
        break;
      }
    }
    if (!titleEl) continue;
    el.setAttribute('title', titleEl.textContent);
    titleEl.remove();
    el.setAttribute('data-bs-toggle', 'tooltip');
    el.setAttribute('data-bs-placement', 'top');
    if (customClass) {
      el.setAttribute('data-bs-custom-class', customClass);
    }
  }
}

// Hides any open tooltip on a promoted node when its modal opens or
// closes. Bootstrap 5.0.2 ties "hide on modal close" to the modal the
// node was inside *at construction* (§3.4) — these nodes start in the
// thumbnail with no modal ancestor, so that built-in hook never fires
// for them, and a move into/out of the modal can otherwise leave a
// tooltip anchored to a stale position. The promoted nodes live in
// either the thumbnail or the modal body (bindSvgModalMoveRestore
// relocates them between the two), never both at once, so both ids
// are checked on every call rather than assuming which one currently
// holds them.
function hideSvgTooltipsAroundModal(thumbId, modalBodyId, itemSelector, modalId) {
  var thumb = document.getElementById(thumbId);
  var modalBody = document.getElementById(modalBodyId);
  var modal = document.getElementById(modalId);
  if (!thumb || !modalBody || !modal) return;
  var hideAll = function() {
    var items = thumb.querySelectorAll(itemSelector);
    var modalItems = modalBody.querySelectorAll(itemSelector);
    for (var i = 0; i < items.length; i++) {
      var instance = bootstrap.Tooltip.getInstance(items[i]);
      if (instance) instance.hide();
    }
    for (var j = 0; j < modalItems.length; j++) {
      var modalInstance = bootstrap.Tooltip.getInstance(modalItems[j]);
      if (modalInstance) modalInstance.hide();
    }
  };
  modal.addEventListener('show.bs.modal', hideAll);
  modal.addEventListener('hide.bs.modal', hideAll);
}
