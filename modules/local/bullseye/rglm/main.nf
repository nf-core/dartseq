process BULLSEYE_R_GLM {
    tag "$meta.id"
    label 'process_medium'

    container params.bullseye_r_container ?: null

    input:
    tuple val(meta), path(cov_tsv), path(mut_tsv), path(coldata_tsv)

    output:
    tuple val(meta), path("${prefix}.glm.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"

    if (params.bullseye_glm_mock) {
        """
        cat > ${prefix}.glm.tsv << 'EOF'
feature	pvalue	lfc
mock_site	1	0
EOF

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye_r: "mock"
        END_VERSIONS
        """
    } else {
        """
        cat > run_bullseye_glm.R << 'EOF'
        suppressPackageStartupMessages({
            library(data.table)
            library(dplyr)
            library(aod)
        })

        source("${params.bullseye_functions_r_file}")

        cov <- fread("${cov_tsv}")
        mut <- fread("${mut_tsv}")
        coldata <- fread("${coldata_tsv}")

        if (!"sample" %in% colnames(coldata)) {
            stop("coldata must contain a 'sample' column")
        }

        rownames(coldata) <- coldata\$sample

        result <- Bullseye(
            data_cov = as.data.frame(cov),
            data_mut = as.data.frame(mut),
            colData = as.data.frame(coldata),
            design = as.formula("${params.bullseye_glm_design}"),
            min.cov = ${params.bullseye_glm_min_cov},
            return.df.only = TRUE,
            link = "${params.bullseye_glm_link}"
        )

        fwrite(as.data.table(result), file = "${prefix}.glm.tsv", sep = "\t")
        EOF

        Rscript run_bullseye_glm.R

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye_r: "local"
        END_VERSIONS
        """
    }
}
