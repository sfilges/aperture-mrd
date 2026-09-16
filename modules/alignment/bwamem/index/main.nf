process BWA_INDEX {
    tag "$fasta"
    label 'process_high_memory'

    container 'community.wave.seqera.io/library/bwa_htslib_samtools:fea458eb782fa887'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
    tuple val("${task.process}"), val('bwa'), eval("bwa 2>&1 | sed -n 's/^Version: //p'"), emit: versions_bwa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''

    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${prefix} \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    // BWA_MEM locates the index by globbing for *.amb, so the stub has to lay down
    // the real file set or a mixed stubbed/real run resolves an empty index path.
    """
    mkdir bwa
    touch bwa/${prefix}.amb
    touch bwa/${prefix}.ann
    touch bwa/${prefix}.bwt
    touch bwa/${prefix}.pac
    touch bwa/${prefix}.sa
    """
}
