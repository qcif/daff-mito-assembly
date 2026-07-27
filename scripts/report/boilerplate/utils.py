"""General utility functions."""

import base64
import re
from pathlib import Path


def serialize(obj):
    """Serialize an object to a JSON string."""
    if hasattr(obj, 'to_json'):
        return obj.to_json()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON"
                    " serializable")


def existing_path(path):
    """Check if a path exists and return a Path object."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path '{path.absolute()}' does not exist.")
    return path.absolute()


def path_safe(dirty: str):
    """Make a string safe for use as a filename."""
    return re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F\s]", "_", dirty)


def get_img_src(path):
    """Return the base64 encoded image source as an HTML img src property."""
    ext = path.suffix[1:]
    if ext.lower() == 'svg':
        return path.read_text()
    return (
        f"data:image/{ext};base64,"
        + base64.b64encode(path.read_bytes()).decode()
    )
