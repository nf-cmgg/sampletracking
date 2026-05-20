#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/sampletracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/sampletracking
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLETRACKING  } from './workflows/sampletracking'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_sampletracking_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_sampletracking_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated file containing information about the samples in the experiment.
    input: Path

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: Path

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // The maximum Y/X reads ratio for a sample to be considered female
    max_yx_female: Float = 0.06

    // The minimum Y/X reads ratio for a sample to be considered male
    min_yx_male: Float = 0.09

    // The minimum coverage of the SRY gene (chrY) for a sample to be considered male
    min_sry_cov_male: Integer = 20

    // The minimum heterozygous SNP fraction for a sample to be considered female
    min_hetx_female: Float = 0.25

    // The maximum heterozygous SNP fraction for a sample to be considered male
    max_hetx_male: Float = 0.05

    // Path to FASTA genome file.
    fasta: Path

    // Path to FASTA index file.
    fai: Path

    // Path to bwa index directory.
    bwa_index: Path

    // Path to haplotype map file.
    haplotype_map: Path

    // Display version and exit.
    version: Boolean = false

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean = false

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean = false

    // Incoming hook URL for messaging service
    hook_url: String = false

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String = 'https://raw.githubusercontent.com/nf-core/test-datasets/'

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String

    // Display the help message.
    help = false

    // Display the full detailed help message.
    help_full: Boolean = false

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean = false
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    SAMPLETRACKING (
        PIPELINE_INITIALISATION.out.samplesheet,
        channel.value([ [id: "bwa"], params.bwa_index ]),
        channel.value([ [id:"genome_fasta"], params.fasta, params.fai ]),
        channel.value([ [id:"haplotype_map"], params.haplotype_map ]),
        params.outdir,
        params.multiqc_config
            ? [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), params.multiqc_config]
            : [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)],
        params.multiqc_logo ? params.multiqc_logo : [],
        params.multiqc_methods_description ? params.multiqc_methods_description : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true),
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        SAMPLETRACKING.out.multiqc_report.filter { meta, _file -> meta.id == "multiqc" }.first(),
    )

    publish:
    multiqc_report = SAMPLETRACKING.out.multiqc_report
    crosscheck_metrics = SAMPLETRACKING.out.crosscheck_metrics
    sex_prediction = SAMPLETRACKING.out.sex_prediction
}

output {
    multiqc_report {
        path { meta, _file ->
            return ("${meta.id}/")
        }
    }
    crosscheck_metrics {
        path { meta, _file ->
            return (meta.pool ? "${meta.pool}/" : "crosscheck_metrics/")
        }
    }
    sex_prediction {
        path { meta, _xy, _sry, _hetx ->
            return (meta.pool ? "${meta.pool}/" : "sex_prediction/")
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
