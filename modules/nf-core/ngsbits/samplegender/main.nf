process NGSBITS_SAMPLEGENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fbf8cfd89c36e9a18a895066bb1da04b93ef585a593b0821ec7037aba6c03474/data':
        'community.wave.seqera.io/library/ngs-bits:2025_12--958625b0e620100a' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_xy.tsv")   , emit: xy_tsv
    tuple val(meta), path("*_hetx.tsv") , emit: hetx_tsv
    tuple val(meta), path("*_sry.tsv")  , emit: sry_tsv
    tuple val("${task.process}"), val('ngsbits'), eval("SampleGender --version  2>&1 | sed 's/SampleGender //'"), topic: versions, emit: versions_ngsbits

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref = fasta ? "-ref ${fasta}" : ""
    """
    SampleGender \\
        -in ${bam} \\
        -method xy \\
        -out ${prefix}_xy.tsv \\
        ${ref} \\
        ${args} \\
    && \\
    SampleGender \\
        -in ${bam} \\
        -method hetx \\
        -out ${prefix}_hetx.tsv \\
        ${ref} \\
        ${args2} \\
    && \\
    SampleGender \\
        -in ${bam} \\
        -method sry \\
        -out ${prefix}_sry.tsv \\
        ${ref} \\
        ${args3}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_xy.tsv
    touch ${prefix}_hetx.tsv
    touch ${prefix}_sry.tsv
    """
}
