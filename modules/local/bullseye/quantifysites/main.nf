process BULLSEYE_QUANTIFY_SITES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(edited_matrix), path(edited_tbi), path(bed)
    val(bullseye_mock)
    val(bullseye_quant_coverage)

    output:
    tuple val(meta), path("${prefix}.quantified.bed"), emit: bed
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    bullseye_version = bullseye_mock ? 'mock' : 'local'
    if (bullseye_mock) {
        """
        cp ${bed} ${prefix}.quantified.bed
        """
    } else {
        """
        perl -S quantify_sites.pl \\
            --EditedMatrix ${edited_matrix} \\
            --EditedMinCoverage ${bullseye_quant_coverage} \\
            --bed ${bed} \\
            --outfile ${prefix}.quantified.bed \\
            --cpu ${task.cpus}
        """
    }
}
