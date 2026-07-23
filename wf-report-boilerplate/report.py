"""Render a workflow report from output files.

This is the entry point invoked by the workflow's report-emitting process.
It loads results from a directory, builds a template context, and writes a
self-contained HTML file (all CSS/JS/images inlined).

Extend `_get_report_context` to add per-analysis result objects.
"""

import json
import logging
import os
import re
from datetime import datetime
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader

from . import config
from .filters.css_hash import css_hash
from .results import Metadata, RunQC
from .utils import get_img_src, serialize

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
config = config.Config()

TEMPLATE_DIR = Path(__file__).parent / 'templates'
STATIC_DIR = Path(__file__).parent / 'static'
EXCLUDE_JS: list[str] = []


def render(
    result_dir: Path,
    samplesheet_file: Path,
    default_params_file: Path,
    params_file: Path,
    versions_file: Path,
    analyst_name: str = None,
    facility: str = None,
):
    """Render an HTML report to the configured output directory."""
    config.load(result_dir)
    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    j2.filters['css_hash'] = css_hash
    j2.filters['html_id'] = lambda s: re.sub(r'[^A-Za-z0-9_-]', '_', str(s))
    template = j2.get_template('index.html')
    context = _get_report_context(
        samplesheet_file,
        default_params_file,
        params_file,
        versions_file,
        analyst_name,
        facility,
    )

    ctx_path = config.result_dir / 'report_context.json'
    with ctx_path.open('w') as f:
        logger.info(f"Writing report context to {ctx_path}")
        json.dump(context, f, indent=2, default=serialize)

    static_files = _get_static_file_contents()
    rendered_html = template.render(**context, **static_files)

    out_path = config.report_path
    with open(out_path, 'w') as f:
        f.write(rendered_html)
    logger.info(f"HTML document written to {out_path}")


def _get_static_file_contents() -> dict:
    """Return static file contents keyed for injection into the template."""
    static_files = {}
    for root, _, files in os.walk(STATIC_DIR):
        root = Path(root)
        if root.name == 'css':
            static_files['css'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in sorted(files)
            ]
        elif root.name == 'js':
            static_files['js'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in sorted(files)
                if not any(p in f for p in EXCLUDE_JS)
            ]
        elif root.name == 'img':
            static_files['img'] = {
                f: get_img_src(root / f)
                for f in sorted(files)
            }
    return {'static': static_files}


def _get_report_context(
    samplesheet_file: Path,
    default_params_file: Path,
    params_file: Path,
    versions_file: Path,
    analyst_name: str,
    facility: str,
) -> dict:
    """Build the context passed to the Jinja template.

    Add per-analysis result objects here as the pipeline stages mature.
    """
    return {
        'title': config.REPORT.TITLE,
        'subtitle_html': config.REPORT.SUBTITLE,
        'sample_id': config.sample_id,
        'analyst_name': analyst_name or '-',
        'facility': facility or '-',
        'versions': _load_yaml(versions_file),
        'start_time': _get_start_time(),
        'end_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'wall_time': _get_walltime(),
        'metadata': _get_metadata(samplesheet_file),
        'parameters': _get_parameters(default_params_file, params_file),
        'run_qc': _get_run_qc(),
        # Add result-object entries here as new stages are implemented.
    }


def _get_start_time():
    if not config.start_time:
        return None
    return config.start_time.strftime("%Y-%m-%d %H:%M:%S")


def _get_walltime():
    """Wall time since workflow start, as zero-padded H/M/S strings."""
    if not config.start_time:
        return None
    seconds = (datetime.now() - config.start_time).total_seconds()
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return {
        'hours': str(int(hours)).zfill(2),
        'minutes': str(int(minutes)).zfill(2),
        'seconds': str(int(seconds)).zfill(2),
    }


def _get_metadata(samplesheet_file: Path) -> 'Metadata | None':
    """Return the row from samples.csv matching the current sample_id.

    Stub: implement parsing against the wf5 samples.csv schema.
    """
    # TODO: parse samplesheet_file, match on sample_id, return Metadata(row)
    return None


def _load_yaml(path: Path) -> dict:
    if not path or not path.exists():
        return {}
    with path.open() as f:
        return yaml.safe_load(f) or {}


def _get_parameters(
    default_params_file: Path,
    params_file: Path,
) -> dict[str, dict[str, str]]:
    """Merge default and user parameter files into a display-ready dict."""
    defaults = _load_yaml(default_params_file)
    user = _load_yaml(params_file)
    return {
        k: {'default': defaults[k], 'user': user.get(k)}
        for k in defaults
    }


def _get_run_qc() -> 'RunQC | None':
    """Return QC metrics for the current sample.

    Stub: implement once QC output format (NanoPlot summary) is defined.
    """
    # TODO: locate QC summary file via config, parse, return RunQC(row)
    return None
