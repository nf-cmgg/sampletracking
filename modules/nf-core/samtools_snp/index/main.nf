process SAMTOOLS_INDEX_SNP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}_snp.bam"),     emit: bam
    tuple val(meta), path("${meta.id}_snp.bam.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp "$input" "${meta.id}_snp.bam"
    samtools index -@ ${task.cpus} "${meta.id}_snp.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """
}
