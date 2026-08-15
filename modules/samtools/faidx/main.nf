process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    container 'biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.{fa,fasta}") , emit: fa , optional: true
    tuple val(meta), path ("*.fai")        , emit: fai, optional: true
    tuple val(meta), path ("*.gzi")        , emit: gzi, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/^samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        $fasta \\
        $args
    """
}