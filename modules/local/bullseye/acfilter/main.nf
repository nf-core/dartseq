process BULLSEYE_ACFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(bed)
    path(fasta)
    val(bullseye_mock)

    output:
    tuple val(meta), path("${prefix}.AC.bed"), emit: bed
    tuple val(meta), path("${prefix}.all_sites.bed"), emit: all_sites_bed
    tuple val(meta), path("${prefix}.all_sites.sequence.bed"), emit: all_sites
    tuple val(meta), path("${prefix}.AC.counts.tsv"), emit: counts
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    bullseye_version = bullseye_mock ? 'mock' : 'local'
    if (bullseye_mock) {
        """
        cp ${bed} ${prefix}.AC.bed
        cp ${bed} ${prefix}.all_sites.bed
        cp ${bed} ${prefix}.all_sites.sequence.bed
        printf "total_sites\\tac_strict_sites\\n0\\t0\\n" > ${prefix}.AC.counts.tsv
        """
    } else {
        """
        perl -F'\\t' -ane 'if (/^\\#/){next;} if(\$F[5] =~ /\\+/){ \$F[1] -= 1; } else { \$F[2] += 1; } print join("\\t", @F)' ${bed} \\
            | bedtools getfasta -bedOut -s -fi ${fasta} -bed - > ${prefix}.2nt.sequence.bed

        perl -F'\\t' -ane 'if(\$F[5] =~ /\\+/){ \$F[1] += 1; } else { \$F[2] -= 1; } print join("\\t", @F)' ${prefix}.2nt.sequence.bed > ${prefix}.all_sites.sequence.bed

        perl -F'\\t' -ane 'if( \$F[-1] =~ m/^AC\$/i ){chomp; pop @F; print join("\\t", @F), "\\n"}' ${prefix}.all_sites.sequence.bed > ${prefix}.AC.bed

        # Same column shape as AC.bed (2nt-context column stripped) but without the
        # motif filter - lets callers propagate every site downstream when the strict
        # canonical-AC-only filter shouldn't be enforced (see skip_bullseye_acfilter).
        perl -F'\\t' -ane 'chomp; pop @F; print join("\\t", @F), "\\n"' ${prefix}.all_sites.sequence.bed > ${prefix}.all_sites.bed

        total_sites=\$(wc -l < ${prefix}.all_sites.sequence.bed)
        ac_strict_sites=\$(wc -l < ${prefix}.AC.bed)
        printf "total_sites\\tac_strict_sites\\n%s\\t%s\\n" "\$total_sites" "\$ac_strict_sites" > ${prefix}.AC.counts.tsv

        rm -f ${prefix}.2nt.sequence.bed
        """
    }
}
