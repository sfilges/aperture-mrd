process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bcftools:1.24--h487d631_1'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.concat.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.concat.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -1 | sed 's/^bcftools //'"), emit: versions_bcftools, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.concat.vcf.gz
    touch ${prefix}.concat.vcf.gz.tbi
    """
}
