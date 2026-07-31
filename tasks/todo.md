# To-do items

- Refactor integration tests and assertions.sh into an nf-test suite
- At the end, do we extract all barcodes or just one? Probably better to
  extract all.
- Check Flye polishing requirement - maybe not necessary any more
- Think about iterative recruitment for tricky/cryptic samples
- Switch flye to --nano-corr if reads are confirmed error-corrected
- If we pass a genome_size hint, that could possibly be a samplesheet column in
  future
- When the `ORGANELLE_MAP` task is drafted, its rendering pass should walk
  both `path1` and `path2` when `plastid_isoforms/` is present, per
  [spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation-ptgaul-derived).
- Examine script config (e.g. hard-coded config in Python scripts like
  bin_target.py) that should perhaps be set as NF params.
