// Renders a small stacked bar per category: a green segment for the
// "subset" quantity with a grey segment stacked above it for the
// remainder up to the "whole", so the full stack reads as the whole.
// Exact values show on hover. A chart whose inputs are all absent
// renders a muted "not available" line instead of a zero-height bar
// (CONSTITUTION principle 7). The chart waits for its container to
// have a width before calling Plotly, because a hidden tab pane sits
// at width 0 and Plotly would otherwise fall back to an oversized
// default.
function renderNestedBarChart(elId, categories, wholeValues, subsetValues, opts) {
  var el = document.getElementById(elId);
  if (!el || typeof Plotly === 'undefined') return;
  opts = opts || {};
  var fmt = opts.valueFormatter || function(v) { return String(v); };
  var hasValue = function(v) { return v !== null && v !== undefined; };

  var haveAny = (wholeValues || []).some(hasValue) ||
    (subsetValues || []).some(hasValue);
  if (!haveAny) {
    el.innerHTML = '<p class="text-muted small">Not available for this sample.</p>';
    return;
  }

  // A width of 0.5 of the category slot makes the gap between
  // adjacent bars equal to the bar width.
  var barWidth = hasValue(opts.barWidthFraction) ? opts.barWidthFraction : 0.5;

  renderWhenSized(el, function() {
    var subset = categories.map(function(_, i) {
      return hasValue((subsetValues || [])[i]) ? subsetValues[i] : null;
    });
    var remainder = categories.map(function(_, i) {
      var whole = (wholeValues || [])[i], part = subset[i];
      if (!hasValue(whole) && !hasValue(part)) return null;
      if (!hasValue(whole)) return 0;
      if (!hasValue(part)) return whole;
      return Math.max(0, whole - part);
    });
    var subsetName = opts.subsetName || 'Subset';
    var remainderName = opts.remainderName || 'Remainder';
    var hoverText = function(values, name) {
      return values.map(function(v) {
        return hasValue(v) ? (name + ': ' + fmt(v)) : '';
      });
    };

    var traces = [
      {
        x: categories, y: subset, type: 'bar', width: barWidth,
        name: subsetName, marker: { color: '#2ca02c' },
        hovertext: hoverText(subset, subsetName), hoverinfo: 'text',
      },
      {
        x: categories, y: remainder, type: 'bar', width: barWidth,
        name: remainderName, marker: { color: '#d9d9d9' },
        hovertext: hoverText(remainder, remainderName), hoverinfo: 'text',
      },
    ];

    var allValues = [];
    (wholeValues || []).forEach(function(v) { if (hasValue(v)) allValues.push(v); });
    (subsetValues || []).forEach(function(v) { if (hasValue(v)) allValues.push(v); });
    var maxValue = allValues.length ? Math.max.apply(null, allValues) : 1;
    var layout = {
      barmode: 'stack',
      yaxis: {
        title: { text: opts.yAxisTitle || '' },
        range: [0, maxValue * 1.15],
      },
      margin: { t: 10, r: 20, l: 50, b: 30 },
      height: 220,
      showlegend: false,
    };
    Plotly.newPlot(el, traces, layout, { displayModeBar: false, responsive: true }).then(function() {
      new ResizeObserver(function() { Plotly.Plots.resize(el); }).observe(el);
    });
  });
}

// Renders one bar per labelled category, all in a single colour. The
// x-tick labels name each bar, so hover shows the value alone. A
// `null` value leaves that bar out rather than drawing a zero-height
// bar, and a chart with no values at all renders a muted "not
// available" line (CONSTITUTION principle 7). `opts.floors` draws
// labelled dashed horizontal lines, e.g. coverage gate floors --
// `floor.label` should include the value itself (e.g. "Fail 10x"),
// since these lines aren't hoverable; the y-axis range always
// includes them. `opts.logY` puts the y-axis on a log scale.
// `opts.yAxisTitle` shows a titled y-axis on a linear scale
// (otherwise the axis is hidden and the value is hover-only).
function renderLabelledBarChart(elId, labels, values, opts) {
  var el = document.getElementById(elId);
  if (!el || typeof Plotly === 'undefined') return;
  opts = opts || {};
  var fmt = opts.valueFormatter || function(v) { return String(v); };
  var hasValue = function(v) { return v !== null && v !== undefined; };

  if (!values.some(hasValue)) {
    el.innerHTML = '<p class="text-muted small">Not available for this sample.</p>';
    return;
  }

  renderWhenSized(el, function() {
    var trace = {
      x: labels, y: values, type: 'bar',
      marker: { color: opts.color || '#2ca02c' },
      hovertext: values.map(function(v) { return hasValue(v) ? fmt(v) : ''; }),
      hoverinfo: 'text',
    };

    var rangeValues = values.filter(hasValue);
    var shapes = [];
    var annotations = [];
    Object.keys(opts.floors || {}).forEach(function(key) {
      var floor = opts.floors[key];
      if (!floor || !hasValue(floor.value)) return;
      rangeValues.push(floor.value);
      shapes.push({
        type: 'line', xref: 'paper', x0: 0, x1: 1,
        y0: floor.value, y1: floor.value,
        line: { color: floor.color, width: 1.5, dash: 'dash' },
      });
      // Plotly places annotations on a log axis in log10 units.
      // floor.label carries its own value (e.g. "Fail 10x") since
      // these dashed lines have no hover of their own.
      annotations.push({
        xref: 'paper', x: 1, xanchor: 'left', yanchor: 'middle',
        y: opts.logY ? Math.log10(floor.value) : floor.value,
        text: floor.label, showarrow: false,
        font: { size: 13, color: floor.color },
      });
    });

    var maxValue = Math.max.apply(null, rangeValues);
    var yaxis;
    if (opts.logY) {
      var positive = rangeValues.filter(function(v) { return v > 0; });
      var minValue = Math.min.apply(null, positive);
      // Log axis ranges are also given in log10 units. The axis stays
      // visible so bar heights are not misread as proportional. dtick
      // 1 keeps ticks to whole powers of ten -- without it Plotly adds
      // unlabelled-looking minor ticks (2, 5) between each decade.
      yaxis = {
        type: 'log', automargin: true, dtick: 1,
        range: [Math.log10(minValue / 2), Math.log10(maxValue * 2)],
        showline: true, linecolor: 'black', linewidth: 1,
        ticks: 'outside', tickcolor: 'black',
      };
    } else if (opts.yAxisTitle) {
      yaxis = {
        title: { text: opts.yAxisTitle }, automargin: true,
        range: [0, maxValue * 1.15],
        showline: true, linecolor: 'black', linewidth: 1,
        ticks: 'outside', tickcolor: 'black',
      };
    } else {
      // Hidden on a linear scale with no axis title: bar heights read
      // on their own, and the exact figure is in the hover text.
      yaxis = { visible: false, range: [0, maxValue * 1.15] };
    }
    var layout = {
      xaxis: {
        title: { text: opts.xAxisTitle || '' }, automargin: true, tickangle: -90,
        showline: true, linecolor: 'black', linewidth: 1,
        ticks: 'outside', tickcolor: 'black', tickfont: { size: 14 },
      },
      yaxis: yaxis,
      margin: { t: 10, r: annotations.length ? 80 : 10, l: 10, b: 10 },
      height: 220,
      showlegend: false,
      shapes: shapes,
      annotations: annotations,
    };
    Plotly.newPlot(el, [trace], layout, { displayModeBar: false, responsive: true }).then(function() {
      new ResizeObserver(function() { Plotly.Plots.resize(el); }).observe(el);
    });
  });
}

// Calls fn(pxWidth) once the element has a non-zero rendered width —
// immediately if it already does, or on the first ResizeObserver
// callback after that (e.g. once its tab pane is shown) otherwise.
function renderWhenSized(el, fn) {
  if (el.clientWidth > 0) {
    fn(el.clientWidth);
    return;
  }
  var observer = new ResizeObserver(function() {
    if (el.clientWidth > 0) {
      observer.disconnect();
      fn(el.clientWidth);
    }
  });
  observer.observe(el);
}
