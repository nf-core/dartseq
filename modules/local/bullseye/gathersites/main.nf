process BULLSEYE_GATHER_SITES {
    tag "$meta.id"
    label 'process_low'

    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(bed_files, stageAs: 'beds??/*')

    output:
    tuple val(meta), path("${prefix}.score.txt"), optional: true, emit: score
    tuple val(meta), path("${prefix}.coverage.txt"), optional: true, emit: coverage
    tuple val(meta), path("${prefix}.mut.txt"), optional: true, emit: mut
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    def outScore = params.bullseye_gather_score ? true : false
    def outCoverage = params.bullseye_gather_coverage ? true : false
    def outMutations = params.bullseye_gather_mutations ? true : false

    if (!outScore && !outCoverage && !outMutations) {
        error("BULLSEYE_GATHER_SITES requires at least one of --bullseye_gather_score, --bullseye_gather_coverage, or --bullseye_gather_mutations")
    }

    if (params.bullseye_mock) {
        """
        header="#chr\tstart\tend\tcluster\tstrand\t${meta.id}"

        if [[ "${outScore}" == "true" ]]; then
            printf "%s\n" "\$header" > ${prefix}.score.txt
            printf "chr1\t100\t101\tmock\t+\t0\n" >> ${prefix}.score.txt
        fi

        if [[ "${outCoverage}" == "true" ]]; then
            printf "%s\n" "\$header" > ${prefix}.coverage.txt
            printf "chr1\t100\t101\tmock\t+\t1\n" >> ${prefix}.coverage.txt
        fi

        if [[ "${outMutations}" == "true" ]]; then
            printf "%s\n" "\$header" > ${prefix}.mut.txt
            printf "chr1\t100\t101\tmock\t+\t1\n" >> ${prefix}.mut.txt
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "mock"
        END_VERSIONS
        """
    } else {
        def gatherFlags = [
            outScore ? '--score' : '',
            outCoverage ? '--coverage' : '',
            outMutations ? '--mutations' : ''
        ].findAll { flag -> flag }.join(' ')

        """
        perl ${params.bullseye_code_dir}/quant/gather_sites.pl \\
            ${gatherFlags} \\
            --outfile ${prefix}.txt \\
            beds*/*.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "local"
        END_VERSIONS
        """
    }
}
