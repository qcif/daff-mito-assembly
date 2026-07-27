// C1 preflight — runs once per invocation before per-sample fan-out.
// Validates the CSV holistically and emits a normalised JSON array.
// A non-zero exit from parse_samplesheet.py aborts the whole run cleanly.

process VALIDATE_SAMPLESHEET {
    label 'process_low'

    input:
    path samplesheet
    path data_dir

    output:
    path 'samples.normalised.json', emit: json

    stub:
    """
    echo '[{"sample_id":"STUB-01","kingdom":"plant","reads":[],"sample_info":"","sample_type":"","sample_receipt_date":"","storage_location":""}]' \
        > samples.normalised.json
    """

    script:
    """
    parse_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        --data-dir    ${data_dir} \\
        --out         samples.normalised.json
    """
}
