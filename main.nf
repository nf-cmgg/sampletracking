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
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCMGG_SAMPLETRACKING {

    take:
    samplesheet // channel: samplesheet read in from --input
    bwa_index   // channel: bwa index file
    genome_fasta // channel: genome fasta file
    haplotype_map // channel: haplotype map file
    outdir       // channel: output directory

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SAMPLETRACKING (
        samplesheet,
        bwa_index,
        genome_fasta,
        haplotype_map,
        outdir
    )
    emit:
    multiqc_report = SAMPLETRACKING.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
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
    NFCMGG_SAMPLETRACKING (
        PIPELINE_INITIALISATION.out.samplesheet,
        channel.value([
            [id: "bwa"],
            file(params.bwa_index, checkIfExists: true)
        ]),
        channel.value(
            [[id:"genome_fasta"],
            file(params.fasta, checkIfExists: true),
            file(params.fai, checkIfExists: true),
        ]),
        channel.value(
            [[id:"haplotype_map"],
            file(params.haplotype_map, checkIfExists: true)
        ]),
        params.outdir
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
        params.hook_url,
        NFCMGG_SAMPLETRACKING.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
