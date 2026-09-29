// Moves an inline-SVG thumbnail's DOM node into its fullscreen modal on
// open and restores it on close, instead of duplicating the markup.
// Two live copies of one inline SVG collide on internal `id`s and
// `url(#...)` references, so the node is relocated, not copied — the
// organelle map (task 44) and the assembly graph (task 47) both need
// this, hence the shared helper.
function bindSvgModalMoveRestore(modalId, thumbId, bodyId) {
  var modal = document.getElementById(modalId);
  var thumb = document.getElementById(thumbId);
  var body = document.getElementById(bodyId);
  if (!modal || !thumb || !body) return;
  modal.addEventListener('show.bs.modal', function() {
    while (thumb.firstChild) body.appendChild(thumb.firstChild);
  });
  modal.addEventListener('hidden.bs.modal', function() {
    while (body.firstChild) thumb.appendChild(body.firstChild);
  });
}
