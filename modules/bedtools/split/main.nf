process BEDTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    tuple val(meta), path(bed)
    val count

    output:
    tuple val(meta), path("*.bed"), emit: beds
    tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed -e "s/bedtools v//g"'), emit: versions_bedtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bedtools split \\
        $args \\
        -n ${count} \\
        -i ${bed} \\
        -p ${prefix}
    """
}