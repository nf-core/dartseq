process BULLSEYE_RACFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(bed)
    path(fasta)
    val(bullseye_mock)
    val(bullseye_rac_remove_seq)

    output:
    tuple val(meta), path("${prefix}.RAC.bed"), emit: bed
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    bullseye_version = bullseye_mock ? 'mock' : 'local'
    if (bullseye_mock) {
        """
        cp ${bed} ${prefix}.RAC.bed
        """
    } else {
        def keepMotifOnly = bullseye_rac_remove_seq ? true : false
        """
        perl -F'\\t' -ane 'if (/^\\#/){next;} if(\$F[5] =~ /\\+/){ \$F[1] -= 2; } else { \$F[2] += 2; } print join("\\t", @F)' ${bed} \\
            | bedtools getfasta -bedOut -s -fi ${fasta} -bed - > ${prefix}.3nt.sequence.bed

        if [[ "${keepMotifOnly}" == "true" ]]; then
            perl -F'\\t' -ane 'if(\$F[5] =~ /\\+/){ \$F[1] += 2; } else { \$F[2] -= 2; } if( \$F[-1] =~ m/[AG]AC/i ){chomp; pop @F; print join("\\t", @F), "\\n"}' ${prefix}.3nt.sequence.bed > ${prefix}.RAC.bed
        else
            perl -F'\\t' -ane 'if(\$F[5] =~ /\\+/){ \$F[1] += 2; } else { \$F[2] -= 2; } if( \$F[-1] =~ m/[AG]AC/i ){ print join "\\t", @F }' ${prefix}.3nt.sequence.bed > ${prefix}.RAC.bed
        fi

        rm -f ${prefix}.3nt.sequence.bed
        """
    }
}
