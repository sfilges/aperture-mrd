process BEDTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed -e "s/bedtools v//g"'), emit: versions_bedtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${beds} > ${prefix}.concat.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.concat.bed
    """
}
