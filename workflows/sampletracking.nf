/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_MEM                       } from '../modules/nf-core/bwa/mem/main'
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
    ch_bwa_index        // channel: [meta, /path/to/bwa_index]
    ch_fasta_fai        // channel: [meta,/path/to/fasta, /path/to/fasta.fai]
    ch_haplotype_map    // channel: [meta, /path/to/haplotype_map]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
    .branch { meta, sample_bam, sample_bam_index, snp_fastq, snp_bam, snp_bam_index ->
        aligned: snp_bam
            return [meta, sample_bam, sample_bam_index, snp_bam, snp_bam_index]
        to_align : snp_fastq
            return [meta, sample_bam, sample_bam_index, snp_fastq]
    }
    .set{ ch_inputs }

    ch_inputs.to_align.multiMap{ meta, sample_bam, sample_bam_index, snp_fastq ->
        fastq:  [meta, snp_fastq]
        bam:    [meta, sample_bam, sample_bam_index]
    }
    .set{ ch_to_align }

    BWA_MEM(
        ch_to_align.fastq,
        ch_bwa_index,
        ch_fasta_fai.map{meta, fasta, fai -> [meta, fasta]},
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    ch_to_align.bam
    .join(BWA_MEM.out.cram, failOnMismatch:true, failOnDuplicate:true)
    .join(BWA_MEM.out.crai, failOnMismatch:true, failOnDuplicate:true)
    .mix(ch_inputs.aligned)
    .map{ meta, sample_bam, sample_bam_index, snp_bam, snp_bam_index ->
        // workaround for bug in MQC crosscheckfingerprints module
        // https://github.com/MultiQC/MultiQC/issues/2449
        return [[id: "crosscheckfingerprints"], sample_bam, sample_bam_index, snp_bam, snp_bam_index]
        //return [[id: meta.pool], sample_bam, sample_bam_index, snp_bam, snp_bam_index]
    }
    .groupTuple()
    .merge(ch_haplotype_map.map{it[1]})
    .map{ meta, sample_bam, sample_bam_index, snp_bam, snp_bam_index, haplotype_map ->
        return [meta, sample_bam.flatten(), sample_bam_index.flatten(), snp_bam.flatten(), snp_bam_index.flatten(), haplotype_map]
    }
    .dump(tag: "Samples to fingerprint", pretty: true)
    .set{ch_to_fingerprint}

    PICARD_CROSSCHECKFINGERPRINTS(
        ch_to_fingerprint,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(PICARD_CROSSCHECKFINGERPRINTS.out.versions)
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
    multiqc_report      = MULTIQC.out.report.toList()                                       // channel: /path/to/multiqc_report.html
    crosscheck_metrics  = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics.map{it[1]}   // channel: [ path(crosscheck_metrics.txt) ]
    versions            = ch_versions                                                       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
