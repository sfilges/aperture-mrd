process PARABRICKS_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1"

    input:
    tuple val(meta), path(input_bam), path(input_index_bam), path(fasta), path(fai)
    
    output:
    tuple val(meta), path("*.metrics.txt") , emit: metrics
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    collectMultipleMetrics \\
        --input $input_bam \\
        --reference $fasta \\
        --output ${prefix}.metrics.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$( collectMultipleMetrics --version )
    END_VERSIONS
    """
}