process FASTQC {
    tag "$meta.id with $task.cpus cores"
    label 'process_medium'

    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    tuple val("${task.process}"), val('fastqc'), eval("fastqc --version | sed '/FastQC v/!d; s/.*v//'"), emit: versions_fastqc, topic: versions

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc \\
        $args \\
        --threads $task.cpus \\
        $reads
    """
}