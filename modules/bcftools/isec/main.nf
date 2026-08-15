process BCFTOOLS_ISEC {
    tag "${meta.id}"
    label 'process_single'

    container 'biocontainers/bcftools:1.24--h487d631_1'

    input:
    tuple val(meta), path(vcfs), path(tbis)
    val consensus_min_count

    output:
    tuple val(meta), path("*.isec.vcf.gz"), emit: vcf
    tuple val(meta), path("*.isec.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -1 | sed 's/^bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_consensus"
    def nfiles = consensus_min_count ? "-n ${consensus_min_count}" : ''
    """
    bcftools isec \\
        ${nfiles} \\
        -w 1 \\
        -Oz -o ${prefix}.isec.vcf.gz \\
        ${args} \\
        ${vcfs.join(' ')}

    tabix -p vcf ${prefix}.isec.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_consensus"
    """
    touch ${prefix}.isec.vcf.gz
    touch ${prefix}.isec.vcf.gz.tbi
    """
}
