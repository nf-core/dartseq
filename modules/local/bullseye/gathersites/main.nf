process BULLSEYE_GATHER_SITES {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(bed_files, stageAs: 'beds??/*')
    val(bullseye_mock)
    val(bullseye_gather_score)
    val(bullseye_gather_coverage)
    val(bullseye_gather_mutations)

    output:
    tuple val(meta), path("${prefix}.score.txt"), optional: true, emit: score
    tuple val(meta), path("${prefix}.coverage.txt"), optional: true, emit: coverage
    tuple val(meta), path("${prefix}.mut.txt"), optional: true, emit: mut
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    bullseye_version = bullseye_mock ? 'mock' : 'local'
    def outScore = bullseye_gather_score ? true : false
    def outCoverage = bullseye_gather_coverage ? true : false
    def outMutations = bullseye_gather_mutations ? true : false

    if (!outScore && !outCoverage && !outMutations) {
        error("BULLSEYE_GATHER_SITES requires at least one of --bullseye_gather_score, --bullseye_gather_coverage, or --bullseye_gather_mutations")
    }

    if (bullseye_mock) {
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
        """
    } else {
        // gather_scores.pl (the vendored Bullseye tool) merges per-sample beds by
        // pushing each sample's value onto an array with NO index tracking - if a
        // site is missing from one sample's bed (e.g. it didn't meet
        // --EditedMinCoverage in QUANTIFY_SITES for that sample, which routinely
        // differs per sample), the row just has fewer columns with no record of
        // WHICH sample dropped out. That's an unrecoverable, ambiguous file
        // downstream (confirmed via a real run: fread() on the result mis-parses
        // the header and run_bullseye_glm.R's rownames<- then fails), and
        // Bullseye() itself (bullseye_functions.R) expects a dense matrix, zero-
        // filling anything below min.cov. So this does the merge as a proper outer
        // join on chr|start|end|strand instead, zero-filling (NA for score, since
        // a missing sample by construction has cov=0 here too) any sample absent
        // at a site - matching what Bullseye() already assumes.
        """
        cat > merge_sites.awk << 'AWK_EOF'
BEGIN {
    FS = "\\t"
    OFS = "\\t"
    nsamples = 0
    for (i = 1; i < ARGC; i++) {
        f = ARGV[i]
        nparts = split(f, parts, "/")
        base = parts[nparts]
        sub(/\\.bed\$/, "", base)
        nsamples++
        names[nsamples] = base
        fileIdx[f] = nsamples
    }
    header = "#chr" OFS "start" OFS "end" OFS "cluster" OFS "strand"
    for (i = 1; i <= nsamples; i++) header = header OFS names[i]
    # Body tmp files are otherwise only created lazily inside the END loop below,
    # which never runs if zero sites are found across all samples - leaving sort
    # with nothing to read. Touch them here so an all-empty result is a valid
    # header-only file instead of a crash.
    if (outCoverage == "true") { print header > covOut; printf "" > covBody }
    if (outMutations == "true") { print header > mutOut; printf "" > mutBody }
    if (outScore == "true") { print header > scoreOut; printf "" > scoreBody }
}
FNR == 1 { curIdx = fileIdx[FILENAME] }
/^#/ { next }
NF < 11 { next }
{
    key = \$1 "|" \$2 "|" \$3 "|" \$6
    if (!(key in seen)) {
        seen[key] = 1
        chrOf[key] = \$1
        startOf[key] = \$2
        endOf[key] = \$3
        clusterOf[key] = \$4
        strandOf[key] = \$6
    }
    covVal[key, curIdx] = \$11
    mutVal[key, curIdx] = \$9
    scoreVal[key, curIdx] = (\$11 == 0 ? "NA" : \$5)
}
END {
    for (key in seen) {
        line = chrOf[key] OFS startOf[key] OFS endOf[key] OFS clusterOf[key] OFS strandOf[key]
        covLine = line
        mutLine = line
        scoreLine = line
        for (i = 1; i <= nsamples; i++) {
            covLine = covLine OFS (((key, i) in covVal) ? covVal[key, i] : 0)
            mutLine = mutLine OFS (((key, i) in mutVal) ? mutVal[key, i] : 0)
            scoreLine = scoreLine OFS (((key, i) in scoreVal) ? scoreVal[key, i] : "NA")
        }
        if (outCoverage == "true") print covLine > covBody
        if (outMutations == "true") print mutLine > mutBody
        if (outScore == "true") print scoreLine > scoreBody
    }
}
AWK_EOF

        awk -v outScore=${outScore} -v outCoverage=${outCoverage} -v outMutations=${outMutations} \\
            -v covOut=${prefix}.coverage.txt -v mutOut=${prefix}.mut.txt -v scoreOut=${prefix}.score.txt \\
            -v covBody=cov.body.tmp -v mutBody=mut.body.tmp -v scoreBody=score.body.tmp \\
            -f merge_sites.awk beds*/*.bed

        if [[ "${outCoverage}" == "true" ]]; then
            LC_ALL=C sort -k1,1 -k2,2n cov.body.tmp >> ${prefix}.coverage.txt
        fi
        if [[ "${outMutations}" == "true" ]]; then
            LC_ALL=C sort -k1,1 -k2,2n mut.body.tmp >> ${prefix}.mut.txt
        fi
        if [[ "${outScore}" == "true" ]]; then
            LC_ALL=C sort -k1,1 -k2,2n score.body.tmp >> ${prefix}.score.txt
        fi

        rm -f merge_sites.awk cov.body.tmp mut.body.tmp score.body.tmp
        """
    }
}
