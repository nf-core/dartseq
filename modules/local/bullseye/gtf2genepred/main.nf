process BULLSEYE_GTF2GENEPRED {
    tag "$meta.id"
    label 'process_medium'

    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(gtf)

    output:
    path "*.refFlat", emit: refFlat
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${gtf.simpleName}"

    """
    perl ${params.bullseye_code_dir}/gtf2genepred.pl \
        --gtf ${gtf} \
        --out ${prefix}.refFlat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bullseye: "local"
    END_VERSIONS
    """
}
