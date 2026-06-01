process BULLSEYE_FIND_EDIT_SITES {
    tag "$meta.id vs $control_meta.id"
    label 'process_medium'

    container params.bullseye_container ?: null

    input:
    tuple val(meta), path(edited_matrix), path(edited_tbi), val(control_meta), path(control_matrix), path(control_tbi), path(annotation_file)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Use contrast-specific parameters if available, otherwise use global params
    def contrast = meta.contrast_params ?: [:]
    def min_edit = contrast.min_edit ?: params.bullseye_min_edit
    def max_edit = contrast.max_edit ?: params.bullseye_max_edit
    def fold_threshold = contrast.fold_threshold ?: params.bullseye_edit_fold_threshold
    def min_sites = contrast.min_sites ?: params.bullseye_min_edit_sites
    def contrast_id = meta.contrast_id ?: 'default'

    prefix = task.ext.prefix ?: "${meta.id}.vs.${control_meta.id}.${contrast_id}"

    if (params.bullseye_mock) {
        """
        cat > ${prefix}.bed << 'EOF'
        chr1\t100\t101\tmock_site
        EOF

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "mock"
        END_VERSIONS
        """
    } else {
        """
        perl ${params.bullseye_code_dir}/Find_edit_site.pl \\
            --annotationFile ${annotation_file} \\
            --EditedMatrix ${edited_matrix} \\
            --controlMatrix ${control_matrix} \\
            --editType ${params.bullseye_edit_type} \\
            --minEdit ${min_edit} \\
            --maxEdit ${max_edit} \\
            --editFoldThreshold ${fold_threshold} \\
            --MinEditSites ${min_sites} \\
            --EditedMinCoverage ${params.bullseye_edited_min_coverage} \\
            --ControlMinCoverage ${params.bullseye_control_min_coverage} \\
            --outfile ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bullseye: "local"
        END_VERSIONS
        """
    }
}
