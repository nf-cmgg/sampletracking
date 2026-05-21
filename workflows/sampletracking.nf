/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_MEM                       } from '../modules/nf-core/bwa/mem/main'
include { PICARD_CROSSCHECKFINGERPRINTS } from '../modules/nf-core/picard/crosscheckfingerprints/main'
include { NGSBITS_SAMPLEGENDER          } from '../modules/nf-core/ngsbits/samplegender/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { samplesheetToList             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_sampletracking_pipeline'
include { SAMTOOLS_INDEX                } from '../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SNP_BAM } from '../modules/nf-core/samtools/index'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLETRACKING {

    take:
    ch_samplesheet              // channel: samplesheet read in from --input
    ch_bwa_index                // channel: [meta, /path/to/bwa_index]
    ch_fasta_fai                // channel: [meta,/path/to/fasta, /path/to/fasta.fai]
    ch_haplotype_map            // channel: [meta, /path/to/haplotype_map]

    outdir                      // path:  path/to/outdir
    multiqc_config              // path:  path/to/multiqc_config
    multiqc_logo                // path:  path/to/multiqc_logo
    multiqc_methods_description // path:  path/to/multiqc_methods_description

    main:
    def ch_multiqc_files = channel.empty()

    def (ch_sample, ch_snp, ch_fastq) = ch_samplesheet
        .multiMap {meta, sample_bam, sample_bam_index, snp_fastq, snp_bam, snp_bam_index ->
            sample: [meta, sample_bam, sample_bam_index]
            snp   : [meta, snp_bam, snp_bam_index]
            fastq : [meta, snp_fastq]
        }

    def (ch_sample_with_idx, ch_sample_no_idx) = ch_sample.branch {
        _meta, _sample_bam, sample_bam_index ->
            with_index: sample_bam_index
            no_index  : !sample_bam_index
    }

    SAMTOOLS_INDEX(
        ch_sample_no_idx.map { meta, sample_bam, _sample_bam_index -> [meta, sample_bam] }
    )

    def ch_sample_fixed = ch_sample_with_idx.mix(
        ch_sample_no_idx
            .join(SAMTOOLS_INDEX.out.index, by: 0)
            .map { meta, sample_bam, _old_index, new_index -> [meta, sample_bam, new_index] }
    )

    def (ch_snp_with_idx, ch_snp_no_idx, ch_snp_none) = ch_snp.branch {
        _meta, snp_bam, snp_bam_index ->
            none      : !snp_bam
            with_index: snp_bam_index
            no_index  : !snp_bam_index
    }

    SAMTOOLS_INDEX_SNP_BAM(
        ch_snp_no_idx.map { meta, snp_bam, _snp_bam_index -> [meta, snp_bam] }
    )

    def ch_snp_fixed = ch_snp_with_idx.mix(
        ch_snp_no_idx
            .join(SAMTOOLS_INDEX_SNP_BAM.out.index, by: 0)
            .map { meta, snp_bam, _old_index, new_index -> [meta, snp_bam, new_index] }
    ).mix(
        ch_snp_none
    )

    def ch_samplesheet_fixed = ch_fastq
        .join(ch_sample_fixed, by: 0)
        .join(ch_snp_fixed,    by: 0)
        .map { meta, snp_fastq, sample_bam, sample_bam_index, snp_bam, snp_bam_index ->
            [meta, sample_bam, sample_bam_index, snp_fastq, snp_bam, snp_bam_index]
        }

    //
    // Crosscheck fingerprints
    //

    def ch_crosscheck_metrics_out = channel.empty()
    ch_samplesheet_fixed
        .filter { meta, _sample_bam, _sample_bam_index, snp_fastq, snp_bam, _snp_bam_index ->
            if(!snp_bam && !snp_fastq) {
                log.warn("No SNP BAM/CRAM/FASTQ files were detected for '${meta.id}'. Skipping the crosscheck fingerprints step for this sample.")
                return false
            }
            return true
        }
        .branch { meta, sample_bam, sample_bam_index, snp_fastq, snp_bam, snp_bam_index ->
            aligned: snp_bam
                return [meta, sample_bam, sample_bam_index, snp_bam, snp_bam_index]
            to_align : snp_fastq
                return [meta, sample_bam, sample_bam_index, snp_fastq]
        }
        .set{ ch_inputs }

    ch_inputs.to_align.multiMap { meta, sample_bam, sample_bam_index, snp_fastq ->
            fastq: [meta, snp_fastq]
            bam:   [meta, sample_bam, sample_bam_index]
        }
        .set { ch_to_align }


    BWA_MEM(
        ch_to_align.fastq,
        ch_bwa_index,
        ch_fasta_fai.map{meta, fasta, _fai -> [meta, fasta]},
        true
    )

    ch_to_align.bam
    .join(BWA_MEM.out.cram, failOnMismatch:true, failOnDuplicate:true)
    .join(BWA_MEM.out.crai, failOnMismatch:true, failOnDuplicate:true)
    .mix(ch_inputs.aligned)
    .map{ meta, sample_bam, sample_bam_index, snp_bam, snp_bam_index ->
        return [groupKey([id: meta.pool], meta.pool_count), sample_bam, sample_bam_index, snp_bam, snp_bam_index]
    }
    .groupTuple()
    .merge(ch_haplotype_map.map{ _meta, haplotype_map -> haplotype_map})
    .map{ meta, sample_bam, sample_bam_index, snp_bam, snp_bam_index, haplotype_map ->
        return [meta, sample_bam.flatten(), sample_bam_index.flatten(), snp_bam.flatten(), snp_bam_index.flatten(), haplotype_map]
    }
    .dump(tag: "Samples to fingerprint", pretty: true)
    .set{ch_to_fingerprint}


    PICARD_CROSSCHECKFINGERPRINTS(
        ch_to_fingerprint,
        ch_fasta_fai
    )
    ch_crosscheck_metrics_out = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics)


    //
    // Determine sample sex
    //
    def ch_sex_prediction_out = channel.empty()
    ch_samplesheet_fixed
        .map { meta, sample_bam, sample_bam_index, _snp_fastq, _snp_bam, _snp_bam_index ->
            [ meta, sample_bam, sample_bam_index ]
        }
        .set { ch_samplegender_input }

    NGSBITS_SAMPLEGENDER(
        ch_samplegender_input,
        ch_fasta_fai
    )

    ch_sex_prediction_out = NGSBITS_SAMPLEGENDER.out.xy_tsv
        .join(NGSBITS_SAMPLEGENDER.out.sry_tsv, failOnMismatch: true, failOnDuplicate: true)
        .join(NGSBITS_SAMPLEGENDER.out.hetx_tsv, failOnMismatch: true, failOnDuplicate: true)

    ch_sex_prediction_out
        .map { meta, xy, sry, hetx ->
            def tsv_list = [:]
            if(!workflow.stubRun) {
                tsv_list = xy.splitCsv(header: ['file', 'sex_yx', 'reads_y', 'reads_x', 'ratio_yx'], sep: "\t", skip:1)[0] +
                    sry.splitCsv(header: ['file', 'sex_sry', 'coverage_sry'], sep: "\t", skip:1)[0] +
                    hetx.splitCsv(header: ['file', 'sex_hetx', 'snps_usable', 'hom_count', 'het_count', 'het_fraction'], sep: "\t", skip:1)[0]
            } else {
                tsv_list = [file:"", sex_xy:"others", reads_y :100, reads_x :100, ratio_xy: 1, sex_sry: 'others', coverage_sry: 10, sex_hetx: 'others', snps_usable: '0 of 0', hom_count: 0, het_count: 0, het_fraction: 0]
            }
            def calculated_sex = "U"
            def predicted_list = [tsv_list.sex_hetx, tsv_list.sex_sry, tsv_list.sex_yx]
            def male_certainty = predicted_list.count("male")/predicted_list.size()
            def female_certainty = predicted_list.count("female")/predicted_list.size()
            def match_certainty = meta.sex == "M" ? male_certainty : meta.sex == "F" ? female_certainty : 0
            if(male_certainty > 0.5) {
                calculated_sex = "M"
            } else if(predicted_list.count("female") >= 2) {
                calculated_sex = "F"
            }
            def data = [
                sample: meta.samplename,
                match_certainty: (match_certainty*100).toInteger(),
                expected_sex: meta.sex,
                calculated_sex: calculated_sex,
                calculated_sex_yx: convert_sex(tsv_list.sex_yx),
                calculated_sex_sry: convert_sex(tsv_list.sex_sry),
                calculated_sex_hetx: convert_sex(tsv_list.sex_hetx),
                ratio_yx: tsv_list.ratio_yx,
                coverage_sry: tsv_list.coverage_sry,
                het_fraction: tsv_list.het_fraction,
            ]
            return [groupKey([id: meta.pool], meta.pool_count), data]
        }
        .groupTuple()
        .map { meta, data ->
            def config = [
                id: "sex_prediction",
                section_name: "Sex prediction",
                description: "This table show the results of the sex prediction with the expected sex and some metrics",
                plot_type: "table",
                headers: [
                    match_certainty: [
                        max: 100,
                        min: 0,
                        suffix: "%"
                    ],
                    expected_sex: [:],
                    calculated_sex: [:],
                    calculated_sex_yx: [:],
                    calculated_sex_sry: [:],
                    calculated_sex_hetx: [:],
                    ratio_yx: [
                        format: "{:,.4f}"
                    ],
                    coverage_sry: [:],
                    het_fraction: [:],
                ],
                data:data.collectEntries { obj -> [obj.sample, obj - obj.subMap("sample")]}
            ]
            return [meta, new org.yaml.snakeyaml.Yaml().dump(config)]
        }
        .collectFile({ meta, content -> ["${meta.id}.sex_prediction_mqc.yml", content] }, newLine:true)
        .map { file ->
            [ [id:file.name.replace(".sex_prediction_mqc.yml", "")], file ]
        }
        .set { ch_sex_prediction_configs }

    ch_multiqc_files = ch_multiqc_files.mix(ch_sex_prediction_configs)

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
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir.toUriString()}/pipeline_info",
            name: 'nf_cmgg_sampletracking_software_mqc_versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    // summary files without meta, e.g. versions, params
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description = channel.value(multiqc_methods_description ? methodsDescriptionText(multiqc_methods_description) : "")

    ch_summary_files = channel.empty()
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        .toList()
        .map { files -> [files] }
        .dump(tag: "Summary files for MultiQC", pretty: true)
    ch_multiqc_input = ch_multiqc_files
        .map { meta, files ->
            def new_meta = meta.pool ? [id: meta.pool] : [id: 'multiqc']
            return [new_meta, files]
        }
        .groupTuple(by: 0)
        .combine(ch_summary_files)
        .map { meta, multiqc_files, summary_files ->
            return [meta, (multiqc_files + summary_files).flatten(), multiqc_config.flatten(), multiqc_logo, [], []]
        }
        .dump(tag: "MULTIQC files", pretty: true)

    MULTIQC (ch_multiqc_input)

    emit:
    multiqc_report      = MULTIQC.out.report        // channel: path(html)
    crosscheck_metrics  = ch_crosscheck_metrics_out // channel: [ val(meta), path(metrics) ]
    sex_prediction      = ch_sex_prediction_out     // channel: [ val(meta), path(tsv) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def convert_sex(input) {
    def result = [
        "male": "M",
        "female": "F"
    ][input]
    return result ?: "U"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
