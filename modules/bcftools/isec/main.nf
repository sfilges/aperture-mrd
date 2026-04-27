process BCFTOOLS_ISEC {
    tag "${meta.id}"
    label 'process_single'

    container 'biocontainers/bcftools:1.21--h8b25389_1'

    publishDir "${params.outdir}/${workflow.runName}/variant_calling/intersection/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcfs), path(tbis)
    val consensus_min_count

    output:
    tuple val(meta), path("*.isec.vcf.gz"), emit: vcf
    tuple val(meta), path("*.isec.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_consensus"
    """
    bcftools isec \\
        -n ${consensus_min_count} \\
        -w 1 \\
        -Oz -o ${prefix}.isec.vcf.gz \\
        ${args} \\
        ${vcfs.join(' ')}

    tabix -p vcf ${prefix}.isec.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_consensus"
    """
    touch ${prefix}.isec.vcf.gz
    touch ${prefix}.isec.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
