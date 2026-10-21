"""General utility functions."""

import base64


def get_img_src(path):
    """Return the base64 encoded image source as an HTML img src property."""
    ext = path.suffix[1:]
    if ext.lower() == 'svg':
        return path.read_text()
    return (
        f"data:image/{ext};base64,"
        + base64.b64encode(path.read_bytes()).decode()
    )


def file_data_uri(path, mime: str) -> str:
    """Return a `data:` URI for an arbitrary file — the self-containment
    mechanism (spec §6a.1) for artefacts that are downloaded or opened
    in a new tab rather than rendered inline (bandage PNG, the
    NanoPlot HTML report, the annotation GFF)."""
    encoded = base64.b64encode(path.read_bytes()).decode()
    return f"data:{mime};base64,{encoded}"


BYTE_UNITS = (('GB', 10 ** 9), ('MB', 10 ** 6), ('KB', 10 ** 3))
SIG_FIGS = 3


def _three_sig_figs(value: float) -> str:
    """Format `value` to three significant figures."""
    digits = SIG_FIGS - len(str(int(value)))
    return f'{value:.{max(digits, 0)}f}'


def format_bytes(n: int) -> str:
    """Return a byte count as decimal (SI) units to three significant
    figures, e.g. `7.54 MB`. The unit is chosen after rounding, so
    999,950,000 reads `1.00 GB`, not `1000 MB`."""
    for i, (unit, factor) in enumerate(BYTE_UNITS):
        if n < factor:
            continue
        text = _three_sig_figs(n / factor)
        if float(text) >= 1000 and i > 0:
            unit, factor = BYTE_UNITS[i - 1]
            text = _three_sig_figs(n / factor)
        return f'{text} {unit}'
    return f'{n} B'
