process R_LIMMA {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(cov_tsv), path(mut_tsv), path(observations_tsv)
    val(bullseye_glm_mock)
    path(bullseye_functions_r_file)
    val(bullseye_glm_design)
    val(bullseye_glm_min_cov)

    output:
    tuple val(meta), path("${prefix}.limma.tsv"), emit: tsv
    tuple val("${task.process}"), val('bullseye_r'), val("${bullseye_r_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye_limma"
    bullseye_r_version = bullseye_glm_mock ? 'mock' : 'local'

    if (bullseye_glm_mock) {
        """
        cat > ${prefix}.glm.tsv << 'EOF'
feature	pvalue	lfc
mock_site	1	0
EOF
        """
    } else {
        """
        cat > run_bullseye_limma.R << 'EOF'
        suppressPackageStartupMessages({
            library(data.table)
            library(dplyr)
            library(limma)
        })

        source("${bullseye_functions_r_file}")

        cov <- fread("${cov_tsv}")
        mut <- fread("${mut_tsv}")
        observations <- as.data.frame(fread("${observations_tsv}"))

        if (!"sample" %in% colnames(observations)) {
            stop("observations must contain a 'sample' column")
        }

        # rownames<- only survives on a data.frame - data.table silently drops it,
        # so observations is converted above before this assignment (Bullseye()
        # requires rownames(colData) to match the sample columns below).
        rownames(observations) <- observations\$sample

        # Bullseye() expects data_cov/data_mut with rownames = site ID and one
        # column per sample (matching rownames(colData)) - GATHER_SITES output
        # instead keeps chr/start/end/cluster/strand as columns, so build the site
        # ID from those and drop them before handing the matrices to Bullseye().
        meta_cols <- c("#chr", "start", "end", "cluster", "strand")
        site_id <- paste(cov[["#chr"]], cov\$start, cov\$end, cov\$strand, sep = "_")

        data_cov <- as.data.frame(cov)[, setdiff(colnames(cov), meta_cols), drop = FALSE]
        data_mut <- as.data.frame(mut)[, setdiff(colnames(mut), meta_cols), drop = FALSE]
        rownames(data_cov) <- site_id
        rownames(data_mut) <- site_id

        # Identify any rows where the coverage is below the specified glm coverage threshold
        rows_w_low_coverage <- rownames(data_cov)[apply(data_cov < ${bullseye_glm_min_cov}, 1, any, na.rm = TRUE)]

        # Remove those rows from both the coverage and mutation dataframes
        data_cov <- data_cov[!rownames(data_cov) %in% rows_w_low_coverage, , drop = FALSE]
        data_mut <- data_mut[!rownames(data_mut) %in% rows_w_low_coverage, , drop = FALSE]
        length(rows_w_low_coverage)   # how many sites were removed?


        # M-values: log2 odds of editing, with a small offset to avoid log(0)
        M <- log2((data_mut + 0.5) / (data_cov - data_mut + 0.5))

        # design for limma
        design <- model.matrix(as.formula("${bullseye_glm_design}"), data = observations)   

        # precision weights from binomial variance of the logit: ~ n * p * (1 - p)
        p <- (data_mut + 0.5) / (data_cov + 1)
        w <- as.matrix(data_cov * p * (1 - p))

        fit <- lmFit(as.matrix(M), design, weights = w)
        fit <- eBayes(fit)
        limmaScores <- topTable(fit, coef = colnames(fit\$coefficients)[2], number = Inf)

        fwrite(as.data.table(limmaScores), file = "${prefix}.limma.tsv", sep = "\t")
        EOF

        Rscript run_bullseye_limma.R
        """
    }
}