//
// Subworkflow with functionality specific to the nf-core/dartseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/dartseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/dartseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    // (absent when starting from matrices via --bullseye_glm_coverage/--bullseye_glm_mut)
    //
    if (!input) {
        ch_samplesheet = channel.empty()
    } else {
        channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map {
                row ->
                    def meta = row[0]
                    def fastq_1 = row.size() > 1 ? row[1] : null
                    def fastq_2 = row.size() > 2 ? row[2] : null
                    def group = row.size() > 3 ? row[3] : null
                    def meta_with_group = group ? meta + [ group: group ] : meta
                    if (!fastq_2) {
                        return [ meta.id, meta_with_group + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta_with_group + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    validatePipelineModeParameters()
}

//
// True when the pipeline starts from coverage/mutation matrices instead of FASTQs
//
def isGlmMatrixMode() {
    return params.bullseye_glm_coverage || params.bullseye_glm_mut
}

def validatePipelineModeParameters() {
    if (isGlmMatrixMode()) {
        if (!params.bullseye_glm_coverage || !params.bullseye_glm_mut) {
            error("Starting from matrices requires both --bullseye_glm_coverage and --bullseye_glm_mut")
        }
        if (!params.bullseye_glm_observations) {
            error("Starting from matrices requires --bullseye_glm_observations")
        }
        if (params.input) {
            error("Use either --input (start from FASTQs) or --bullseye_glm_coverage/--bullseye_glm_mut (start from matrices), not both")
        }
        if (!params.bullseye_glm_mock && !params.bullseye_functions_r_file) {
            error("Starting from matrices requires --bullseye_functions_r_file unless --bullseye_glm_mock is enabled")
        }
        // FASTQ, alignment and Bullseye site-calling checks below do not apply
        return
    }

    if (!params.input) {
        error("--input is required unless starting from matrices (--bullseye_glm_coverage, --bullseye_glm_mut, --bullseye_glm_observations)")
    }

    def allowed_trimmers = ['none', 'fastp', 'trimgalore']
    if (!allowed_trimmers.contains(params.trimmer as String)) {
        error("Invalid --trimmer '${params.trimmer}'. Allowed values: ${allowed_trimmers.join(', ')}")
    }

    def allowed_aligners = ['star', 'hisat2']
    if (!allowed_aligners.contains(params.aligner as String)) {
        error("Invalid --aligner '${params.aligner}'. Allowed values: ${allowed_aligners.join(', ')}")
    }

    if (!params.skip_alignment && params.aligner == 'star') {
        if (!params.star_index && !params.fasta) {
            error("STAR alignment selected but neither --star_index nor --fasta is set (or derivable from --genome)")
        }
        if (!params.gtf) {
            error("STAR alignment selected but --gtf is not set (or not derivable from --genome)")
        }
    }

    if (!params.skip_alignment && params.aligner == 'hisat2') {
        if (!params.hisat2_index && !params.fasta) {
            error("HISAT2 alignment selected but neither --hisat2_index nor --fasta is set")
        }
        if (!params.hisat2_splicesites && !params.gtf) {
            error("HISAT2 alignment selected but neither --hisat2_splicesites nor --gtf is set")
        }
    }

    if (!params.skip_bullseye && params.skip_alignment) {
        error("--skip_bullseye false requires alignment outputs. Please disable --skip_alignment.")
    }

    if (!params.skip_bullseye && !params.bullseye_mock && !params.annotation_reFlat_file && !params.gtf) {
        error("--skip_bullseye false requires either --annotation_reFlat_file or --gtf unless --bullseye_mock is enabled")
    }

    if (!params.skip_bullseye_racfilter && params.skip_bullseye) {
        error("--skip_bullseye_racfilter false requires --skip_bullseye false")
    }

    // BULLSEYE_ACFILTER (mandatory A-before-edited-C check) always runs whenever Bullseye
    // does, so --fasta is required unconditionally now - not just for the optional RAC filter.
    if (!params.skip_bullseye && !params.bullseye_mock && !params.fasta) {
        error("--skip_bullseye false requires --fasta unless --bullseye_mock is enabled")
    }

    if (!params.skip_bullseye_gather_sites && params.skip_bullseye) {
        error("--skip_bullseye_gather_sites false requires --skip_bullseye false")
    }

    if (!params.skip_bullseye_gather_sites && !(params.bullseye_gather_score || params.bullseye_gather_coverage || params.bullseye_gather_mutations)) {
        error("--skip_bullseye_gather_sites false requires at least one of --bullseye_gather_score, --bullseye_gather_coverage, or --bullseye_gather_mutations")
    }

    if (!params.skip_bullseye_glm && params.skip_bullseye_gather_sites) {
        error("--skip_bullseye_glm false requires --skip_bullseye_gather_sites false")
    }

    if (!params.skip_bullseye_glm && !params.bullseye_glm_mock && !params.bullseye_functions_r_file) {
        error("--skip_bullseye_glm false requires --bullseye_functions_r_file unless --bullseye_glm_mock is enabled")
    }

    if (!params.skip_bullseye_glm && !params.bullseye_glm_observations) {
        error("--skip_bullseye_glm false requires --bullseye_glm_observations")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Add in-text citations for tools used by this pipeline.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Function is rendered in methodsDescriptionText for MultiQC.
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Bibliographic entries for tools used by this pipeline.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Function is rendered in methodsDescriptionText for MultiQC.
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
