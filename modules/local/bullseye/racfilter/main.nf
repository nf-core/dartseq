process BULLSEYE_RACFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(bed)
    path(fasta)

    output:
    tuple val(meta), path("${prefix}.RAC.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    if (params.bullseye_mock) {
        """
        cp ${bed} ${prefix}.RAC.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "mock"
        END_VERSIONS
        """
    } else {
        def keepMotifOnly = params.bullseye_rac_remove_seq ? true : false
        """
        perl -F'\\t' -ane 'if (/^\\#/){next;} if(\$F[5] =~ /\\+/){ \$F[1] -= 2; } else { \$F[2] += 2; } print join("\\t", @F)' ${bed} \\
            | bedtools getfasta -bedOut -s -fi ${fasta} -bed - > ${prefix}.3nt.sequence.bed

        if [[ "${keepMotifOnly}" == "true" ]]; then
            perl -F'\\t' -ane 'if(\$F[5] =~ /\\+/){ \$F[1] += 2; } else { \$F[2] -= 2; } if( \$F[-1] =~ m/[AG]AC/i ){chomp; pop @F; print join("\\t", @F), "\\n"}' ${prefix}.3nt.sequence.bed > ${prefix}.RAC.bed
        else
            perl -F'\\t' -ane 'if(\$F[5] =~ /\\+/){ \$F[1] += 2; } else { \$F[2] -= 2; } if( \$F[-1] =~ m/[AG]AC/i ){ print join "\\t", @F }' ${prefix}.3nt.sequence.bed > ${prefix}.RAC.bed
        fi

        rm -f ${prefix}.3nt.sequence.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "local"
        END_VERSIONS
        """
    }
}
