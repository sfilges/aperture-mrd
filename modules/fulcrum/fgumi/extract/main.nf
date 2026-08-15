process FULCRUM_FGUMI_EXTRACT{
    // https://github.com/fulcrumgenomics/fgumi
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/fgumi:0.4.0--h54198d6_0'

    input:
    tuple val(meta), path(reads)
    path(fasta) // Reference FASTA file (must be indexed with .fai). Required. Required when error is selected
    val(library)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi extract \\
        --inputs ${reads.join(' ')} \\
        --output ${prefix}.bam \\
        ${args} \\
        --sample ${prefix} \\
        --library "${library}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}