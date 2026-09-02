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
