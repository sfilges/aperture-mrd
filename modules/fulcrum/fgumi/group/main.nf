process FULCRUM_FGUMI_GROUP {
    tag {meta.id}
    label 'process_medium'

    input:
    tuple val(meta), path(sorted_bam)
    val(umi_mode)

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // 'paired' for duplex or 'adjacency' for simplex/codec workflows
    def grouping_strategy = umi_mode == 'duplex' ? 'paired' : 'adjacency'
    """
    fgumi group \\
        --input ${sorted_bam} \\
        --output ${prefix}_grouped.bam  \\
        --strategy ${grouping_strategy} \\
        --threads ${task.cpus} \\
        ${args}
    """
}