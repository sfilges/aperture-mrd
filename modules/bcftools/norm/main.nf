process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bcftools:1.24--h487d631_1'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.norm.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.norm.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -1 | sed 's/^bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools norm \\
        -m -both \\
        -f $fasta \\
        $args \\
        -Oz -o ${prefix}.norm.vcf.gz \\
        $vcf

    tabix -p vcf ${prefix}.norm.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.norm.vcf.gz
    touch ${prefix}.norm.vcf.gz.tbi
    """
}
