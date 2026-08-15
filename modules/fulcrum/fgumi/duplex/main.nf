process FULCRUM_FGUMI_DUPLEX {
    tag {meta.id}
    label 'process_high'

    container 'quay.io/biocontainers/fgumi:0.4.0--h54198d6_0' 

    input:
    tuple val(meta), path(grouped_bam)
    val min_reads
    val keep_rejected

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.stats.txt"), emit: stats, optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rejects_command = keep_rejected ? "--rejects ${prefix}.rejects.bam" : ''
    """
    fgumi duplex \\
        --input ${grouped_bam} \\
        --output ${prefix}.bam \\
        --min-reads ${min_reads} \\
        --threads ${task.cpus} \\
        --stats ${prefix}.stats.txt \\
        ${rejects_command} \\
        ${args}
    """
}