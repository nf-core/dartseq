process BULLSEYE_QUANTIFY_SITES {
    tag "$meta.id"
    label 'process_medium'

    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(edited_matrix), path(edited_tbi), path(bed)

    output:
    tuple val(meta), path("${prefix}.quantified.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    if (params.bullseye_mock) {
        """
        cp ${bed} ${prefix}.quantified.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "mock"
        END_VERSIONS
        """
    } else {
        """
        perl ${params.bullseye_code_dir}/quant/quantify_sites.pl \\
            --EditedMatrix ${edited_matrix} \\
            --EditedMinCoverage ${params.bullseye_quant_coverage} \\
            --bed ${bed} \\
            --outfile ${prefix}.quantified.bed \\
            --cpu ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "local"
        END_VERSIONS
        """
    }
}
