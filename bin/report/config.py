"""Static report constants.

Task 42 §5.4 inlines every diagnostic payload into `metadata.json`, so
the renderer's result *values* come from that single JSON document
(see `report.py`), never from globbing a per-sample result directory.
This module carries no result-file accessors — the prior glob-based
`Config` class (`_get_file_by_pattern`, `sample_id`, `flags`,
`default_params`) pointed at files that never existed in this repo and
has been retired (task 43a §3).
"""

REPO_URL = 'https://github.com/qcif/daff-biosecurity-wf5'
WF_NAME = 'WF5 (Organelle Assembly + Barcode Recovery)'

REPORT_TITLE = "Organelle assembly + barcode recovery report"
REPORT_SUBTITLE_HTML = (
    "Results generated through the"
    f' <a href="{REPO_URL}" target="_blank">'
    f'{WF_NAME} Nextflow pipeline</a>.'
)
