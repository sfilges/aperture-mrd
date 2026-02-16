process BCFTOOLS_ISEC {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bcftools:1.12.0--h8b12537_0'

    input:
    tuple val(meta), path(bcf_a), path(bcf_b), path(bcf_c)
    val operator
    val n_files

    output:
    tuple val(meta), path('*.bed'), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_intersect"
    """
    bcftools \\
        isec \\
        -p ${prefix} \\
        -n ${operator}${n_files} \\
        -w1 -O z -p \\
        $bcf_a \\
        $bcf_b \\
        $bcf_c \\
        $args \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}