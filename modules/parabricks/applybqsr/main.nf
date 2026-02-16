process PARABRICKS_APPLYBQSR {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1"

    input:
    tuple val(meta), path(input_bam), path(input_index_bam), path(input_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bqsr.bam")        , emit: bqsr_bam
    tuple val(meta), path("*.bqsr.bai")        , emit: bqsr_bai
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parabricks ApplyBQSR \\
        --bam $input_bam \\
        --bai $input_index_bam \\
        --reference $fasta \\
        --output ${prefix}.bqsr.bam \\
        --outputBai ${prefix}.bqsr.bai \\
        $args
     
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$( parabricks --version )
    END_VERSIONS
    """
}