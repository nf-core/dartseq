process RUSTQC {
    tag "$meta.id"
    label 'process_low'
    stageInMode 'copy'
    container { params.rustqc_container ?: null }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.rustqc.tsv"), emit: tsv
    tuple val(meta), path("*.rustqc.json"), emit: json
    tuple val("${task.process}"), val('rustqc'), val(params.rustqc_mock ? 'mock' : 'external'), emit: versions, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (params.rustqc_mock) {
        return """
        printf "metric\tvalue\nmock_reads\t1\n" > ${prefix}.rustqc.tsv
        printf "{\"sample\":\"${prefix}\",\"status\":\"mock\"}\n" > ${prefix}.rustqc.json
        """
    }

    if (!params.rustqc_cmd) {
        error("--run_rustqc requires --rustqc_cmd unless --rustqc_mock is enabled")
    }

    def cmd_template = params.rustqc_cmd as String
    def cmd = (cmd_template.contains('{bam}') || cmd_template.contains('{prefix}'))
        ? cmd_template.replace('{bam}', bam.toString()).replace('{prefix}', prefix)
        : "${cmd_template} ${bam} ${prefix}"

    """
    ${cmd}

    if [[ ! -s ${prefix}.rustqc.tsv ]]; then
        echo "RustQC command did not produce ${prefix}.rustqc.tsv" >&2
        exit 1
    fi

    if [[ ! -s ${prefix}.rustqc.json ]]; then
        printf '{"sample":"%s","status":"generated_from_tsv"}\n' "${prefix}" > ${prefix}.rustqc.json
    fi
    """
}
