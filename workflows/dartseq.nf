/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { HISAT2_ALIGN           } from '../modules/nf-core/hisat2/align/main'
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD           } from '../modules/nf-core/hisat2/build/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BULLSEYE_PARSEBAM      } from '../modules/local/bullseye/parsebam/main'
include { BULLSEYE_GTF2GENEPRED  } from '../modules/local/bullseye/gtf2genepred/main'
include { BULLSEYE_FIND_EDIT_SITES } from '../modules/local/bullseye/findeditsites/main'
include { BULLSEYE_SUMMARIZE_SITES } from '../modules/local/bullseye/summarizesites/main'
include { BULLSEYE_QUANTIFY_SITES } from '../modules/local/bullseye/quantifysites/main'
include { BULLSEYE_ACFILTER      } from '../modules/local/bullseye/acfilter/main'
include { BULLSEYE_RACFILTER     } from '../modules/local/bullseye/racfilter/main'
include { BULLSEYE_GATHER_SITES  } from '../modules/local/bullseye/gathersites/main'
include { BULLSEYE_R_GLM         } from '../modules/local/bullseye/rglm/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dartseq_pipeline'
include { isGlmMatrixMode        } from '../subworkflows/local/utils_nfcore_dartseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DARTSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def ch_trimmed_reads = ch_samplesheet
    def ch_aligned_bam = channel.empty()
    def ch_gtf = channel.empty()
    // [ meta, coverage, mut ] matrices for BULLSEYE_R_GLM: from BULLSEYE_GATHER_SITES
    // when starting from FASTQs, or from --bullseye_glm_coverage/--bullseye_glm_mut
    def ch_bullseye_glm_matrices = channel.empty()
    // Starting from matrices: ch_samplesheet is empty, so FastQC/trimming run no
    // tasks, and alignment and Bullseye site calling are skipped below.
    def glm_matrix_mode = isGlmMatrixMode()
    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // Optional read trimming (fastp|trimgalore|none)
    //
    if (params.trimmer == 'fastp') {
        ch_fastp_input = ch_samplesheet.map { meta, reads -> [ meta, reads, [] ] }
        FASTP (
            ch_fastp_input,
            false,
            false,
            false
        )
        ch_trimmed_reads = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.log.collect{it[1]})
    } else if (params.trimmer == 'trimgalore') {
        TRIMGALORE (
            ch_samplesheet
        )
        ch_trimmed_reads = TRIMGALORE.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.html.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})
    }

    //
    // Optional alignment (star|hisat2) → sort → index
    //
    ch_sorted_bam = channel.empty()

    if (!params.skip_alignment && !glm_matrix_mode) {
        ch_fasta = params.fasta ? Channel.value([ [ id: 'fasta' ], file(params.fasta, checkIfExists: true) ]) : channel.empty()
        ch_gtf = params.gtf ? Channel.value([ [ id: 'gtf' ], file(params.gtf, checkIfExists: true) ]) : channel.empty()

        if (params.aligner == 'star') {
            if (!params.star_index && !params.fasta) {
                error "STAR alignment requires either --star_index or --fasta (plus --gtf)."
            }
            if (!params.gtf) {
                error "STAR alignment requires --gtf."
            }

            ch_star_index = params.star_index
                ? Channel.value([ [ id: 'star_index' ], file(params.star_index, checkIfExists: true) ])
                : null

            if (!params.star_index) {
                STAR_GENOMEGENERATE (
                    ch_fasta,
                    ch_gtf
                )
                ch_star_index = STAR_GENOMEGENERATE.out.index
            }

            STAR_ALIGN (
                ch_trimmed_reads,
                ch_star_index,
                ch_gtf,
                params.star_ignore_sjdbgtf
            )
            ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})
            // conf/modules.config pins STAR_ALIGN's ext.args to
            // `--outSAMtype BAM SortedByCoordinate`, so bam_sorted_aligned is the
            // only channel that's ever actually populated for this pipeline. Do
            // NOT .mix() it with .bam/.bam_sorted as a "fallback chain": the
            // module's `bam` emit glob ('*d.out.bam') also matches
            // "*.Aligned.sortedByCoord.out.bam", so the same file gets captured on
            // both channels - mixing them double-counted every sample, cascading
            // into a 2x/4x/16x task-count multiplication all the way down through
            // Bullseye's edited x control pairing.
            ch_sorted_bam = STAR_ALIGN.out.bam_sorted_aligned

        } else if (params.aligner == 'hisat2') {
            if (!params.hisat2_index && !params.fasta) {
                error "HISAT2 alignment requires either --hisat2_index or --fasta (plus --gtf)."
            }

            ch_hisat2_index = params.hisat2_index
                ? Channel.value([ [ id: 'hisat2_index' ], file(params.hisat2_index, checkIfExists: true) ])
                : null

            ch_hisat2_splicesites = params.hisat2_splicesites
                ? Channel.value([ [ id: 'splicesites' ], file(params.hisat2_splicesites, checkIfExists: true) ])
                : null

            if (!params.hisat2_splicesites) {
                if (!params.gtf) {
                    error "HISAT2 without --hisat2_splicesites requires --gtf."
                }
                HISAT2_EXTRACTSPLICESITES (
                    ch_gtf
                )
                ch_hisat2_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt
            }

            if (!params.hisat2_index) {
                if (!params.gtf) {
                    error "HISAT2 index build requires --gtf."
                }
                HISAT2_BUILD (
                    ch_fasta,
                    ch_gtf,
                    ch_hisat2_splicesites
                )
                ch_hisat2_index = HISAT2_BUILD.out.index
            }

            HISAT2_ALIGN (
                ch_trimmed_reads,
                ch_hisat2_index,
                ch_hisat2_splicesites,
                false
            )
            ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]})

            // HISAT2 outputs unsorted BAM → sort it
            SAMTOOLS_SORT (
                HISAT2_ALIGN.out.bam,
                Channel.value([ [ id: 'no_fasta' ], [], [] ]),
                ''
            )
            ch_sorted_bam = SAMTOOLS_SORT.out.bam
        }

        // Index sorted BAMs (required by Bullseye)
        SAMTOOLS_INDEX ( ch_sorted_bam )
        ch_aligned_bam = ch_sorted_bam.join( SAMTOOLS_INDEX.out.index )
    }

    //
    // Optional Bullseye steps
    //
    if (!params.skip_bullseye && !glm_matrix_mode) {
        // Stage annotation into Bullseye tasks so it is always visible inside containers.
        ch_bullseye_annotation = channel.empty()
        if (params.bullseye_mock) {
            ch_bullseye_annotation = channel.value(file("$projectDir/assets/bullseye.mock.refFlat", checkIfExists: true))
        } else if (params.annotation_reFlat_file) {
            ch_bullseye_annotation = channel.value(file(params.annotation_reFlat_file, checkIfExists: true))
        } else {
            BULLSEYE_GTF2GENEPRED ( ch_gtf )
            ch_bullseye_annotation = BULLSEYE_GTF2GENEPRED.out.refFlat
            // versions now flow via topic: versions (see BULLSEYE_GTF2GENEPRED main.nf) -
            // no manual mix needed, and .out.versions is now a tuple, not a Path.
        }

        BULLSEYE_PARSEBAM ( ch_aligned_bam, params.bullseye_mock, params.bullseye_parse_min_coverage )

        // Two modes: contrast-based (explicit comparisons) or group-based (all edited vs all control)
        if (params.bullseye_contrasts) {
            // Parse contrasts CSV and create specific pairs
            ch_contrasts = channel
                .fromPath(params.bullseye_contrasts, checkIfExists: true)
                .splitCsv(header: true)
                .map { row ->
                    def mode = row.mode ?: 'standard'
                    def min_edit = row.min_edit ?: (mode == 'differential' ? '3' : '5')
                    def max_edit = row.max_edit ?: (mode == 'differential' ? '95' : '90')
                    def fold_threshold = row.fold_threshold ?: (mode == 'differential' ? '1.2' : '1.5')
                    def min_sites = row.min_sites ?: '3'
                    [
                        contrast_id: row.contrast_id,
                        edited_group: row.edited_group,
                        control_group: row.control_group,
                        mode: mode,
                        min_edit: min_edit,
                        max_edit: max_edit,
                        fold_threshold: fold_threshold,
                        min_sites: min_sites
                    ]
                }

            // For each contrast, find matching samples and create pairs
            ch_bullseye_pairs = ch_contrasts
                .combine(BULLSEYE_PARSEBAM.out.matrix)
                .branch {
                    contrast, meta, matrix, tbi ->
                        edited: (meta.group ?: '').toString() == contrast.edited_group
                            return [ contrast, meta, matrix, tbi ]
                        control: (meta.group ?: '').toString() == contrast.control_group
                            return [ contrast, meta, matrix, tbi ]
                        other: true
                }

            // Group edited and control samples per contrast
            ch_bullseye_edited_per_contrast = ch_bullseye_pairs.edited
                .map { contrast, meta, matrix, tbi ->
                    [ contrast.contrast_id, contrast, meta, matrix, tbi ]
                }
                .groupTuple()

            ch_bullseye_control_per_contrast = ch_bullseye_pairs.control
                .map { contrast, meta, matrix, tbi ->
                    [ contrast.contrast_id, contrast, meta, matrix, tbi ]
                }
                .groupTuple()

            // Join and create all pairs within each contrast
            ch_bullseye_pairs = ch_bullseye_edited_per_contrast
                .join(ch_bullseye_control_per_contrast)
                .flatMap { _contrast_id, contrast_list, edited_meta_list, edited_matrix_list, edited_tbi_list,
                           _contrast_list2, control_meta_list, control_matrix_list, control_tbi_list ->
                    def contrast = contrast_list[0]  // All entries have same contrast info
                    def pairs = []
                    edited_meta_list.eachWithIndex { edited_meta, ei ->
                        control_meta_list.eachWithIndex { control_meta, ci ->
                            // Enrich meta with contrast info
                            def enriched_meta = edited_meta + [
                                contrast_id: contrast.contrast_id,
                                contrast_mode: contrast.mode
                            ]
                            pairs << [
                                enriched_meta,
                                edited_matrix_list[ei],
                                edited_tbi_list[ei],
                                control_meta,
                                control_matrix_list[ci],
                                control_tbi_list[ci],
                                contrast
                            ]
                        }
                    }
                    return pairs
                }

            // Extract unique edited samples for quantification (collect all edited_group values from contrasts)
            ch_bullseye_edited = ch_contrasts
                .map { contrast -> contrast.edited_group }
                .unique()
                .combine(BULLSEYE_PARSEBAM.out.matrix)
                .filter { edited_group, meta, _matrix, _tbi ->
                    (meta.group ?: '').toString() == edited_group
                }
                .map { _edited_group, meta, matrix, tbi ->
                    [ meta, matrix, tbi ]
                }

        } else {
            // Legacy mode: all edited vs all control based on group label
            def control_group = (params.bullseye_control_group ?: 'control').toString().toLowerCase()
            ch_bullseye_control = BULLSEYE_PARSEBAM.out.matrix.filter { meta, _matrix, _tbi ->
                (meta.group ?: '').toString().toLowerCase() == control_group
            }
            ch_bullseye_edited = BULLSEYE_PARSEBAM.out.matrix.filter { meta, _matrix, _tbi ->
                (meta.group ?: '').toString().toLowerCase() != control_group
            }

            // Use default contrast parameters
            def default_contrast = [
                contrast_id: 'default',
                mode: 'standard',
                min_edit: params.bullseye_min_edit.toString(),
                max_edit: params.bullseye_max_edit.toString(),
                fold_threshold: params.bullseye_edit_fold_threshold.toString(),
                min_sites: params.bullseye_min_edit_sites.toString()
            ]

            ch_bullseye_pairs = ch_bullseye_edited
                .combine(ch_bullseye_control)
                .map { edited_meta, edited_matrix, edited_tbi, control_meta, control_matrix, control_tbi ->
                    [ edited_meta, edited_matrix, edited_tbi, control_meta, control_matrix, control_tbi, default_contrast ]
                }
        }

        // Add annotation file to all pairs
        ch_bullseye_pairs_with_annotation = ch_bullseye_pairs
            .combine(ch_bullseye_annotation)
            .map { edited_meta, edited_matrix, edited_tbi, control_meta, control_matrix, control_tbi, contrast, annotation_file ->
                [ edited_meta + [ contrast_params: contrast ], edited_matrix, edited_tbi, control_meta, control_matrix, control_tbi, annotation_file ]
            }

        BULLSEYE_FIND_EDIT_SITES (
            ch_bullseye_pairs_with_annotation,
            params.bullseye_mock,
            params.bullseye_min_edit,
            params.bullseye_max_edit,
            params.bullseye_edit_fold_threshold,
            params.bullseye_min_edit_sites,
            params.bullseye_edit_type,
            params.bullseye_edited_min_coverage,
            params.bullseye_control_min_coverage
        )

        ch_bullseye_grouped_sites = BULLSEYE_FIND_EDIT_SITES.out.bed
            .map { meta, bed -> [ meta.id, meta, bed ] }
            .groupTuple()
            .map { _sample_id, metas, beds -> [ metas[0], beds ] }

        BULLSEYE_SUMMARIZE_SITES ( ch_bullseye_grouped_sites, params.bullseye_mock, params.bullseye_replicates_min, params.bullseye_min_mut )

        ch_bullseye_edited_keyed = ch_bullseye_edited
            .map { meta, matrix, tbi -> [ meta.id, meta, matrix, tbi ] }

        ch_bullseye_summarized_keyed = BULLSEYE_SUMMARIZE_SITES.out.bed
            .map { meta, bed -> [ meta.id, meta, bed ] }

        ch_bullseye_quantify_input = ch_bullseye_edited_keyed
            .join(ch_bullseye_summarized_keyed)
            .map { _sample_id, edited_meta, matrix, tbi, _summary_meta, bed ->
                [ edited_meta, matrix, tbi, bed ]
            }

        BULLSEYE_QUANTIFY_SITES ( ch_bullseye_quantify_input, params.bullseye_mock, params.bullseye_quant_coverage )

        ch_bullseye_sites_for_post = BULLSEYE_QUANTIFY_SITES.out.bed

        // Always run ACFILTER: it produces AC.counts.tsv/all_sites*.bed regardless
        // (useful diagnostics either way), but whether the strict AC-only subset is
        // what actually propagates downstream is controlled by skip_bullseye_acfilter -
        // it's not yet confirmed whether APOBEC1-YTH only edits in that canonical
        // A-before-C context, so callers can opt to keep non-canonical sites too.
        ch_bullseye_fasta = channel.value(file(params.fasta, checkIfExists: true))
        BULLSEYE_ACFILTER ( ch_bullseye_sites_for_post, ch_bullseye_fasta, params.bullseye_mock )
        ch_bullseye_sites_for_post = params.skip_bullseye_acfilter
            ? BULLSEYE_ACFILTER.out.all_sites_bed
            : BULLSEYE_ACFILTER.out.bed

        if (!params.skip_bullseye_racfilter) {
            BULLSEYE_RACFILTER ( ch_bullseye_sites_for_post, ch_bullseye_fasta, params.bullseye_mock, params.bullseye_rac_remove_seq )
            ch_bullseye_sites_for_post = BULLSEYE_RACFILTER.out.bed
        }

        if (!params.skip_bullseye_gather_sites) {
            ch_bullseye_gather_input = ch_bullseye_sites_for_post
                .map { _meta, bed -> [ [ id: 'all_bullseye_sites' ], bed ] }
                .groupTuple()
                // groupTuple() collects each sample's bed in whatever order its
                // upstream task happens to finish - not guaranteed stable across
                // runs. BULLSEYE_GATHER_SITES numbers samples by this list's order
                // (stageAs 'beds??/*'), so an unsorted list makes the merged
                // coverage/mut/score matrices' sample-column order (and thus their
                // content hash) non-deterministic, even though the data is the same.
                .map { meta, beds -> [ meta, beds.sort { it.name } ] }

            BULLSEYE_GATHER_SITES (
                ch_bullseye_gather_input,
                params.bullseye_mock,
                params.bullseye_gather_score,
                params.bullseye_gather_coverage,
                params.bullseye_gather_mutations
            )

            if (!params.skip_bullseye_glm) {
                ch_bullseye_glm_matrices = BULLSEYE_GATHER_SITES.out.coverage
                    .join(BULLSEYE_GATHER_SITES.out.mut)
            }
        }
    }

    //
    // Optional Bullseye GLM (differential editing) - on gathered matrices, or
    // directly on user-provided matrices (--bullseye_glm_coverage/--bullseye_glm_mut)
    //
    if (glm_matrix_mode) {
        ch_bullseye_glm_matrices = channel.value([
            [ id: 'glm_matrix' ],
            file(params.bullseye_glm_coverage, checkIfExists: true),
            file(params.bullseye_glm_mut, checkIfExists: true)
        ])
    }

    if (glm_matrix_mode || !params.skip_bullseye_glm) {
        ch_bullseye_glm_observations = channel.value(file(params.bullseye_glm_observations, checkIfExists: true))
        ch_bullseye_glm_functions_r_file = channel.value(file(params.bullseye_functions_r_file, checkIfExists: true))

        ch_bullseye_glm_input = ch_bullseye_glm_matrices
            .combine(ch_bullseye_glm_observations)
            .map { meta, cov, mut, observations -> [ meta, cov, mut, observations ] }

        BULLSEYE_R_GLM (
            ch_bullseye_glm_input,
            params.bullseye_glm_mock,
            ch_bullseye_glm_functions_r_file,
            params.bullseye_glm_design,
            params.bullseye_glm_min_cov,
            params.bullseye_glm_link
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'dartseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'dartseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
