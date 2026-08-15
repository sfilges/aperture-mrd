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
    // Truncate to chrom/start/end: blacklists mix 3- and 4-column BEDs, and a ragged
    // column count makes the downstream `bedtools sort` abort. The name column is
    // unused here — BEDTOOLS_MERGE collapses these to plain intervals anyway.
    """
    awk 'BEGIN { FS = OFS = "\\t" } !/^(#|track|browser)/ && NF >= 3 { print \$1, \$2, \$3 }' ${beds} \\
        > ${prefix}.concat.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.concat.bed
    """
}
