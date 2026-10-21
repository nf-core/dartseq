process BULLSEYE_FIND_EDIT_SITES {
    tag "$meta.id vs $control_meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bullseye:1.0.0--3e5fd1a65e0e10e9'

    input:
    tuple val(meta), path(edited_matrix), path(edited_tbi), val(control_meta), path(control_matrix), path(control_tbi), path(annotation_file)
    val(bullseye_mock)
    val(bullseye_min_edit)
    val(bullseye_max_edit)
    val(bullseye_edit_fold_threshold)
    val(bullseye_min_edit_sites)
    val(bullseye_edit_type)
    val(bullseye_edited_min_coverage)
    val(bullseye_control_min_coverage)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('bullseye'), val("${bullseye_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Use contrast-specific parameters if available, otherwise use global params
    def contrast = meta.contrast_params ?: [:]
    def min_edit = contrast.min_edit ?: bullseye_min_edit
    def max_edit = contrast.max_edit ?: bullseye_max_edit
    def fold_threshold = contrast.fold_threshold ?: bullseye_edit_fold_threshold
    def min_sites = contrast.min_sites ?: bullseye_min_edit_sites
    def contrast_id = meta.contrast_id ?: 'default'

    prefix = task.ext.prefix ?: "${meta.id}.vs.${control_meta.id}.${contrast_id}"
    bullseye_version = bullseye_mock ? 'mock' : 'local'

    if (bullseye_mock) {
        """
        cat > ${prefix}.bed << 'EOF'
        chr1\t100\t101\tmock_site
        EOF
        """
    } else {
        """
        perl -S Find_edit_sites.pl \\
            --annotationFile ${annotation_file} \\
            --EditedMatrix ${edited_matrix} \\
            --controlMatrix ${control_matrix} \\
            --editType ${bullseye_edit_type} \\
            --minEdit ${min_edit} \\
            --maxEdit ${max_edit} \\
            --editFoldThreshold ${fold_threshold} \\
            --MinEditSites ${min_sites} \\
            --EditedMinCoverage ${bullseye_edited_min_coverage} \\
            --ControlMinCoverage ${bullseye_control_min_coverage} \\
            --outfile ${prefix}.bed
        """
    }
}
