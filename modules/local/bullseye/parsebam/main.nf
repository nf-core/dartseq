process BULLSEYE_PARSEBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(bam), path(bai)
    val(bullseye_mock)
    val(bullseye_parse_min_coverage)

    output:
    tuple val(meta), path("${prefix}.matrix.gz"), path("${prefix}.matrix.gz.tbi"), emit: matrix
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    bullseye_version = bullseye_mock ? 'mock' : 'local'
    def removeDuplicates = task.ext.removeDuplicates ?: false
    if (bullseye_mock) {
        """
        cat > ${prefix}.matrix.gz << 'EOF'
        chr1\t100\t+\t10\t2\tA\tG
        EOF

        touch ${prefix}.matrix.gz.tbi
        """
    } else {
        """
        perl -S parseBAM.pl \\
            --input ${bam} \\
            --output ${prefix}.matrix \\
            --cpu ${task.cpus} \\
            --minCoverage ${bullseye_parse_min_coverage} \\
            ${removeDuplicates ? '--removeDuplicates' : ''}

        # parseBAM may emit either <prefix>.matrix or <prefix>.matrix.gz depending on mode/version.
        if [[ -s ${prefix}.matrix ]]; then
            bgzip -f ${prefix}.matrix
        elif [[ ! -s ${prefix}.matrix.gz ]]; then
            echo "parseBAM did not produce ${prefix}.matrix or ${prefix}.matrix.gz" >&2
            exit 1
        fi

        if [[ ! -s ${prefix}.matrix.gz.tbi ]]; then
            tabix -f -s 1 -b 2 -e 2 ${prefix}.matrix.gz
        fi
        """
    }
}
