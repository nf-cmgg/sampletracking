/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM2_MEM                   } from '../modules/nf-core/bwamem2/mem/main'
include { PICARD_CROSSCHECKFINGERPRINTS } from '../modules/nf-core/picard/crosscheckfingerprints/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap              } from 'plugin/nf-validation'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_sampletracking_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLETRACKING {

    take:
    ch_samplesheet      // channel: samplesheet read in from --input
    ch_bwamem2_index    // channel: [meta, /path/to/bwamem2_index]
    ch_fasta            // channel: [meta,/path/to/fasta]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
    .branch { meta, sample_bam, snp_fastq, snp_bam, haplotype_map ->
        aligned: snp_bam
            return [meta, sample_bam, snp_bam, haplotype_map]
        to_align : snp_fastq
            return [meta, sample_bam, snp_fastq, haplotype_map]
    }
    .set{ ch_inputs }

    ch_inputs.to_align.multiMap{ meta, sample_bam, snp_fastq, haplotype_map ->
        fastq:          [meta, snp_fastq]
        haplotype_map:  [meta, sample_bam, haplotype_map]
    }
    .set{ ch_to_align }

    BWAMEM2_MEM(
        ch_to_align.fastq,
        ch_bwamem2_index,
        true
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    ch_to_align.haplotype_map
    .join(BWAMEM2_MEM.out.bam, failOnMismatch:true, failOnDuplicate:true)
    .dump(tag: "algined_snp")
    .map{ meta, sample_bam, haplotype_map, snp_bam ->
        [meta, sample_bam, snp_bam, haplotype_map]}
    .mix(ch_samplesheet.aligned)
    .set(ch_to_fingerprint)

    PICARD_CROSSCHECKFINGERPRINTS(
        ch_to_fingerprint,
        ch_fasta
    )
    ch_versions = ch_versions.mix(PICARD_CROSSCHECKFINGERPRINTS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics.map{it[1]})

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
