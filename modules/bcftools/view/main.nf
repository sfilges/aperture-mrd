process BCFTOOLS_VIEW {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bcftools:1.21--h8b25389_1'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(targets)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_cmd = targets ? "-T ^${targets}" : ''
    """
    bcftools view \\
        ${targets_cmd} \\
        $args \\
        -Oz -o ${prefix}.filtered.vcf.gz \\
        $vcf

    tabix -p vcf ${prefix}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
