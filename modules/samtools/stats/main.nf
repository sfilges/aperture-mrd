process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir "${params.outdir}/${workflow.runName}/reports/samtools/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//'"), emit: versions_samtools, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input} \\
        ${args} \\
        > ${prefix}.stats
    """
}