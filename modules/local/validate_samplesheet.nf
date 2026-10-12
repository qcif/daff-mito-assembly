// C1 preflight — runs once per invocation before per-sample fan-out.
// Validates the CSV holistically and emits a normalised JSON array.
// A non-zero exit from parse_samplesheet.py aborts the whole run cleanly.

process VALIDATE_SAMPLESHEET {
    label 'process_low'

    input:
    path samplesheet
    val data_dir_set

    output:
    path 'samples.normalised.json', emit: json

    script:
    """
    parse_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        ${data_dir_set ? '--data-dir-set' : ''} \\
        --out         samples.normalised.json
    """

    stub:
    """
    echo '[{"sample_id":"STUB-01","assembly_target":"animal_mt","reads":[],"sample_info":"","sample_type":"","sample_receipt_date":"","storage_location":""}]' \
        > samples.normalised.json
    """
}
