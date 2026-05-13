process BULLSEYE_PLACEHOLDER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('bullseye_placeholder'), val('0.1.0'), emit: versions, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Placeholder for Bullseye integration.
    # Replace this process with parseBAM/find_edit_site/summarize_sites/quantify_sites calls.
    printf "chr1\t1\t2\t${prefix}|mock_site\n" > ${prefix}.bullseye_sites.bed
    printf "site\tmetric\n${prefix}_site\t1\n" > ${prefix}.bullseye_summary.tsv
    """
}
