process BWAMEM3_INDEX {
    // https://bwa-mem3.readthedocs.io/en/v0.6.0/
    tag "${meta.id}"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa-mem3_samtools_htslib:c28a809633c294ed'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: aligned
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index, optional: true
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem3 \\
        index \\
        ${fasta}
    """
}