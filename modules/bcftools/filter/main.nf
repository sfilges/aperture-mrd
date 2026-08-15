process BCFTOOLS_FILTER {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bcftools:1.24--h487d631_1'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -1 | sed 's/^bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools filter \\
        $args \\
        -Oz -o ${prefix}.filtered.vcf.gz \\
        $vcf

    tabix -p vcf ${prefix}.filtered.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi
    """
}
