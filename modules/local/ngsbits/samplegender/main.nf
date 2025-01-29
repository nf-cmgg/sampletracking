process NGSBITS_SAMPLEGENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngs-bits:2024_11--py312hd80e9a6_0':
        'biocontainers/ngs-bits:2024_11--py312hd80e9a6_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_xy.tsv")   , emit: xy_tsv
    tuple val(meta), path("*_hetx.tsv") , emit: hetx_tsv
    tuple val(meta), path("*_sry.tsv")  , emit: sry_tsv
    path "versions.yml"                 , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(SampleGender --version 2>&1) | sed 's/SampleGender //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_xy.tsv
    touch ${prefix}_hetx.tsv
    touch ${prefix}_sry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(SampleGender --version 2>&1) | sed 's/SampleGender //' )
    END_VERSIONS
    """
}
