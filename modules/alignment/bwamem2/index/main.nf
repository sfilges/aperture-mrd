process BWAMEM2_INDEX {
    tag "$fasta"
    // NOTE Requires 28N GB memory where N is the size of the reference sequence
    // source: https://github.com/bwa-mem2/bwa-mem2/issues/9
    memory { 28.B * fasta.size() }

    container 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2"), emit: index
    tuple val("${task.process}"), val('bwamem2'), eval("echo \$(bwa-mem2 version 2>&1) | sed 's/.* //'"), emit: versions_bwamem2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir bwamem2
    bwa-mem2 \\
        index \\
        $args \\
        -p bwamem2/${prefix} \\
        $fasta
    """
}