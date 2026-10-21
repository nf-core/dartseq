process BULLSEYE_GTF2GENEPRED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(gtf)

    output:
    path "*.refFlat", emit: refFlat
    tuple val("${task.process}"), val('bullseye'), val('local'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${gtf.simpleName}"

    """
    perl -S gtf2genepred.pl \
        --gtf ${gtf} \
        --out ${prefix}.refFlat
    """
}
