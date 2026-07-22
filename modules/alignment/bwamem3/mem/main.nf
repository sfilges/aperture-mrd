process BWAMEM3_MEM {
    // https://bwa-mem3.readthedocs.io/en/v0.6.0/
    tag "${meta.id}"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa-mem3_htslib_samtools:391ed2ac52c4a15a'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: cram
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index, optional: true
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem3 \\
        mem \\
        -R "${meta.read_group}" \\
        ${args} \\
        -t ${task.cpus} \\
        \$INDEX \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} ${fasta} -O cram -o ${prefix}.cram -
    """
}