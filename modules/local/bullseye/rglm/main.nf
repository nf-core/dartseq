process BULLSEYE_R_GLM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/r-aod_r-base_r-data.table_r-dplyr_pruned:50e53b1dbf9130ec'

    input:
    tuple val(meta), path(cov_tsv), path(mut_tsv), path(observations_tsv)
    val(bullseye_glm_mock)
    path(bullseye_functions_r_file)
    val(bullseye_glm_design)
    val(bullseye_glm_min_cov)
    val(bullseye_glm_link)

    output:
    tuple val(meta), path("${prefix}.glm.tsv"), emit: tsv
    tuple val("${task.process}"), val('bullseye_r'), val("${bullseye_r_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
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
        # aod::betabin's matrix ops (vcov/qr.solve) run through BLAS/LAPACK; a
        # multi-threaded BLAS is not bit-reproducible across thread counts, which
        # otherwise makes the fitted p-values/log2FC flip in their last few digits
        # between runs on the same input.
        export OPENBLAS_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        export MKL_NUM_THREADS=1

        cat > run_bullseye_glm.R << 'EOF'
        suppressPackageStartupMessages({
            library(data.table)
            library(dplyr)
            library(tidyr)   # pivot_longer/nest/replace_na used by Bullseye()
            library(tibble)  # rownames_to_column used by Bullseye()
            library(aod)
        })

        source("${bullseye_functions_r_file}")

        cov <- fread("${cov_tsv}")
        mut <- fread("${mut_tsv}")
        observations <- as.data.frame(fread("${observations_tsv}"))
        design <- as.formula("${bullseye_glm_design}")
        group_col <- all.vars(design)[1]

        if (!"sample" %in% colnames(observations)) {
            stop("observations must contain a 'sample' column")
        }
        if (!group_col %in% colnames(observations)) {
            stop(sprintf("observations must contain the design column '%s'", group_col))
        }

        # Matrices: GATHER_SITES format - chr/start/end/cluster/strand, then one
        # column per sample, identical columns and rows in coverage and mutation.
        meta_cols <- c("#chr", "start", "end", "cluster", "strand")
        if (!all(meta_cols %in% colnames(cov)) || !identical(colnames(cov), colnames(mut))) {
            stop("coverage and mutation matrices must have identical columns: ",
                 paste(meta_cols, collapse = ", "), ", then one column per sample")
        }
        site_id <- paste(cov[["#chr"]], cov\$start, cov\$end, cov\$strand, sep = "_")
        if (!identical(site_id, paste(mut[["#chr"]], mut\$start, mut\$end, mut\$strand, sep = "_"))) {
            stop("coverage and mutation matrices must list the same sites in the same order")
        }

        missing_samples <- setdiff(observations\$sample, setdiff(colnames(cov), meta_cols))
        if (length(missing_samples) > 0) {
            stop("observations samples not found in the matrix columns: ", paste(missing_samples, collapse = ", "))
        }

        groups <- as.character(observations[[group_col]])
        if (length(unique(groups)) < 2) {
            stop(sprintf("design column '%s' needs at least 2 groups to compare, found: %s",
                         group_col, paste(unique(groups), collapse = ", ")))
        }

        # Bullseye() assigns sample columns to groups by regex (".*<group>.*") on the
        # column names, so name each column <group>__<sample> and make sure no group
        # label also matches another group's columns.
        new_names <- paste0(groups, "__", observations\$sample)
        for (lev in unique(groups)) {
            hit <- grepl(paste0(".*", lev, ".*"), new_names, perl = TRUE)
            if (!identical(hit, groups == lev)) {
                stop(sprintf("group '%s' also matches samples of another group (%s) - group labels must not occur inside other group or sample names",
                             lev, paste(new_names[hit & groups != lev], collapse = ", ")))
            }
        }

        if (nrow(cov) == 0) {
            warning("coverage/mutation matrices contain no sites - writing an empty GLM result")
            fwrite(data.table(site = character()), file = "${prefix}.glm.tsv", sep = "\\t")
            quit(save = "no")
        }

        # Bullseye() expects rownames = site ID and only the sample columns, in the
        # same order as rownames(colData). rownames<- only survives on a data.frame
        # (data.table silently drops it), hence the as.data.frame() conversions.
        data_cov <- as.data.frame(cov)[, observations\$sample, drop = FALSE]
        data_mut <- as.data.frame(mut)[, observations\$sample, drop = FALSE]
        colnames(data_cov) <- new_names
        colnames(data_mut) <- new_names
        rownames(data_cov) <- site_id
        rownames(data_mut) <- site_id
        # Bullseye() re-adds the rownames as a 'sample' column, so colData must not
        # carry one already.
        col_data <- observations[, setdiff(colnames(observations), "sample"), drop = FALSE]
        rownames(col_data) <- new_names

        result <- Bullseye(
            data_cov = data_cov,
            data_mut = data_mut,
            colData = col_data,
            design = design,
            min.cov = ${bullseye_glm_min_cov},
            return.df.only = TRUE,
            link = "${bullseye_glm_link}"
        )

        fwrite(as.data.table(result), file = "${prefix}.glm.tsv", sep = "\t")
        EOF

        Rscript run_bullseye_glm.R
        """
    }
}
