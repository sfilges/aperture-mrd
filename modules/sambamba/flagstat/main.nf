process SAMBAMBA_FLAGSTAT {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('sambamba'), eval("echo \$(sambamba --version 2>&1) | awk '{print \$2}'"), emit: versions_sambamba, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sambamba \\
        flagstat \\
        -t $task.cpus \\
        $bam \\
        > ${prefix}.stats
    """
}