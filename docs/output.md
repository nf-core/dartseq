# nf-core/dartseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

This section summarizes the key output folders produced by default profiling modes.
When optional modules are enabled, additional directories such as `bullseye/`
are also written with sample-level result files.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and detect RNA editing sites, using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [Bullseye](#bullseye) - RNA editing site detection and motif filtering (optional)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Bullseye

<details markdown="1">
<summary>Output files</summary>

- `bullseye/<sample>/`
  - `*.bullseye.AC.bed`: Edited sites that pass the A-before-edited-C check (strict `AC` 2-nt motif). This is the set used for downstream steps unless `--skip_bullseye_acfilter` is set.
  - `*.bullseye.all_sites.bed`: Every candidate edited site, same column shape as `AC.bed` but without the motif filter. Used downstream instead of `AC.bed` when `--skip_bullseye_acfilter` is set (e.g. to keep non-canonical sites, since it's not yet confirmed whether APOBEC1-YTH edits only in the canonical `AC` context).
  - `*.bullseye.all_sites.sequence.bed`: Same sites as `all_sites.bed`, annotated with its 2-nt sequence context column. Use this file to look into which non-canonical motifs are present (e.g. sites where the C is upstream of the A, "CA").
  - `*.bullseye.AC.counts.tsv`: Two columns, `total_sites` and `ac_strict_sites` - the total number of candidate conversions detected and how many of those fall in the strict `AC` motif, so the fraction removed by the strict filter is visible per sample regardless of whether it's enforced downstream.
  - `*.bullseye.RAC.bed`: present unless `--skip_bullseye_racfilter` is set. Sites further restricted to the `[A/G]AC` 3-nt motif.
- `bullseye/all_bullseye_sites/`: present unless `--skip_bullseye_gather_sites` is set. Per-sample site-level BED files merged into matrices across all samples.
  - `*.bullseye.coverage.txt`: present when `--bullseye_gather_coverage` is set (default). Coverage at each site, one column per sample; sites missing from a sample are filled with `0`.
  - `*.bullseye.mut.txt`: present when `--bullseye_gather_mutations` is set (default). Mutation (edit) counts at each site, one column per sample; same zero-filling as above.
  - `*.bullseye.score.txt`: present when `--bullseye_gather_score` is set (default). Per-sample editing score at each site; sites with zero coverage in a sample are `NA` rather than `0`.
  - `*.bullseye.glm.tsv`: present unless `--skip_bullseye_glm` is set. Per-site GLM test results (`feature`, `pvalue`, `lfc`) comparing groups defined in `--bullseye_glm_observations`, fit from the `coverage.txt`/`mut.txt` matrices above. See the [Bullseye GLM usage section](usage.md#bullseye-glm) for how to configure this step; `NA` p-values at a site are an expected result of the fitted model (e.g. no residual degrees of freedom), not a pipeline error.

</details>

Bullseye's edited-site detection is followed by a check requiring the base immediately before an edited C to be A (the canonical APOBEC1-YTH editing signal). By default this is enforced (`AC.bed` propagates downstream); set `--skip_bullseye_acfilter` to instead propagate every candidate site (`all_sites.bed`) if there's reason to suspect editing outside that canonical context. The filter step itself always runs regardless, so `AC.counts.tsv` and the two `all_sites*` files are always available to compare against.

Site-level BEDs from every sample can optionally be merged into matrices (`--skip_bullseye_gather_sites false`; skipped by default) and, from those, tested for differential editing with a GLM (`--skip_bullseye_glm false`; also skipped by default).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
