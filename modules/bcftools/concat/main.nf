process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bcftools:1.21--h8b25389_1'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.concat.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.concat.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Index inputs for overlap-aware concat
    for vcf in $vcfs; do
        bcftools index --tbi \$vcf
    done

    bcftools concat \\
        -a \\
        $args \\
        $vcfs \\
        -Ou | bcftools sort -Oz -o ${prefix}.concat.vcf.gz

    tabix -p vcf ${prefix}.concat.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.concat.vcf.gz
    touch ${prefix}.concat.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
