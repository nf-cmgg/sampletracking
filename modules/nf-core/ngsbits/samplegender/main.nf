process NGSBITS_SAMPLEGENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db759890fb18613dd6178305e20a588bda85a12c3d06f885899aca2f54725985/data'
:         'community.wave.seqera.io/library/ngs-bits:2026_06--10de2f01af4c9c32' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*_xy.tsv"), emit: xy_tsv
    tuple val(meta), path("*_hetx.tsv"), emit: hetx_tsv
    tuple val(meta), path("*_sry.tsv"), emit: sry_tsv
    tuple val("${task.process}"), val('ngsbits'), eval("SampleGender --version  2>&1 | sed -n 's/SampleGender //p'"), topic: versions, emit: versions_ngsbits

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
