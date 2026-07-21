process SAMTOOLS_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/samtools:1.19.2--h50ea8bc_0'

    publishDir "${params.outdir}/${workflow.runName}/preprocessing/markduplicates/${meta.id}/", mode: params.publish_dir_mode, pattern: "*.{cram,crai}"
    publishDir "${params.outdir}/${workflow.runName}/reports/", mode: params.publish_dir_mode, pattern: "*.metrics"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.cram"),    emit: cram
    tuple val(meta), path("*.crai"),    emit: crai
    tuple val(meta), path("*.metrics"), emit: metrics
    tuple val("${task.process}"), val('samtools'), eval("echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools collate -O -u --threads ${task.cpus} ${bam} | \\
        samtools fixmate -m -u --threads ${task.cpus} - - | \\
        samtools sort -u --threads ${task.cpus} - | \\
        samtools markdup \\
            $args \\
            -f ${prefix}.md.metrics \\
            --threads ${task.cpus} \\
            --reference ${fasta} \\
            --output-fmt cram \\
            - ${prefix}.sorted.md.cram

    samtools index -@${task.cpus} ${prefix}.sorted.md.cram
    """
}
