process BULLSEYE_SUMMARIZE_SITES {
    tag "$meta.id"
    label 'process_low'

    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(bed_files, stageAs: 'repOnly??/*')

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye.sites"
    if (params.bullseye_mock) {
        """
        cp ${bed_files[0]} ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "mock"
        END_VERSIONS
        """
    } else {
        def repOnlyArgs = bed_files.collect { bedFile -> "--repOnly ${bedFile}" }.join(' ')
        """
        perl ${params.bullseye_code_dir}/summarize_sites.pl \\
            --MinRep ${params.bullseye_replicates_min} \\
            ${repOnlyArgs} > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "local"
        END_VERSIONS
        """
    }
}
