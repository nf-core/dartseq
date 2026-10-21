process BULLSEYE_SUMMARIZE_SITES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(bed_files, stageAs: 'repOnly??/*')
    val(bullseye_mock)
    val(bullseye_replicates_min)
    val(bullseye_min_mut)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye.sites"
    bullseye_version = bullseye_mock ? 'mock' : 'local'
    if (bullseye_mock) {
        """
        cp ${bed_files[0]} ${prefix}.bed
        """
    } else {
        // --repOnly is a boolean switch (merge the same edited sample against several
        // controls without computing coverage/average-editing-rate stats), not a
        // per-file option - the bed files themselves are positional arguments.
        def bedFileArgs = bed_files.collect { bedFile -> "${bedFile}" }.join(' ')
        """
        perl -S summarize_sites.pl \\
            --MinRep ${bullseye_replicates_min} \\
            --mut ${bullseye_min_mut} \\
            --repOnly \\
            ${bedFileArgs} > ${prefix}.bed
        """
    }
}
