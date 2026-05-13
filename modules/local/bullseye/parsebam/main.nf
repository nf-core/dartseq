process BULLSEYE_PARSEBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.22.1 bioconda::htslib=1.22.1"
    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.matrix.gz"), path("${prefix}.matrix.gz.tbi"), emit: matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def removeDuplicates = task.ext.removeDuplicates ?: false
    if (params.bullseye_mock) {
        """
        cat > ${prefix}.matrix.gz << 'EOF'
        chr1\t100\t+\t10\t2\tA\tG
        EOF

        touch ${prefix}.matrix.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "mock"
        END_VERSIONS
        """
    } else {
        """
        perl ${params.bullseye_code_dir}/parseBAM.pl \\
            --input ${bam} \\
            --output ${prefix}.matrix \\
            --cpu ${task.cpus} \\
            --minCoverage ${params.bullseye_parse_min_coverage} \\
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "local"
        END_VERSIONS
        """
    }
}
