process SAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bam")  , emit: bam ,   optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.bai")  , emit: bai ,   optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/^samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"

    """
    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        $args \\
        $input \\
        -o ${prefix}.${output_extension}

    samtools index -@${task.cpus} ${prefix}.${output_extension}
    """
}