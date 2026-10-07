# Improvements to report UI

> Reviewed report: reference-material/task-54-report-int-animal-01.html

## Overview

- "Key findings" should use smaller font for "Top BLAST hit" value
- "Key findings" and "Warnings" should be shown in two columns, 8-4 width respectively
- "Key findings" should replace row colour with colours badges for row values, showing either green tick, orange warn (!), red danger (X) icons in the badge, and place the badge to the left of the text.

### Sequencing data quality

- Should display table and plot in two columns, equal width.
- Table "Clean reads" column title should drop " (post CHOPPER + FILTLONG)" as too verbose
- Plot should switch to vertical bars to save horizontal space


## Validation

### Recruitment and coverage gate

- Move "Gate decision" to it's own large badge component at the top
- Replace some table values with vertical bar charts:
    - Stacked bar chart for aligned/recruited reads in grey/green
    - Stacked bar chart for raw/recruited bases in grey/green
    - Stacked bar chart for est/recruited coverage in grey/green, with coverage floors annotated as horizontal lines
- Move subsampled (last 4 rows) to dedicated table hidden in a modal, activated by a "Subsampling" button, only shown when subsampling has occurred.
- Genetic codes should not show raw ID (e.g. "5") - please map these to human-readable names e.g. "invertebrate"


## Assembly

### Assembly

- Add two buttons for "Assembly FASTA": "copy to clipboard" and "download" (use icons for both)
- The bandage image should be shown in the "Assembly statistics" sections as it has nothing to do with annotation. Please create a third column for this (`col-sm-3`), next to "Per-contig mean coverage".
- The bandage image seems to be displayed really large now, and is overflowing its div and taking over the whole page.
- The "Genome map" title should be removed now, and its info badge moved to the "Genome annotation" section title.

**Genome map**
- Some text at the bottom of the SVG (below the legend) is being cut off. Add some height?
- Would be nice if we could label the barcode loci only - perhaps some text placed outside the circle?


## Barcodes

Don't show "Panel barcodes not recovered" section when all were recovered.
