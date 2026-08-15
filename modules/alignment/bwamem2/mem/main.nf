process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: cram
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index, optional: true
    tuple val("${task.process}"), val('bwamem2'), eval("bwa-mem2 version 2>&1 | tail -n1 | sed 's/.* //'"), emit: versions_bwamem2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/^samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = { params.split_fastq > 1 ? "${meta.id}".concat('.').concat(reads.get(0).name.tokenize('.')[0]) : "${meta.id}.sorted" }()
    def args = task.ext.args ?: ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem2 \\
        mem \\
        -R "${meta.read_group}" \\
        -t $task.cpus \\
        $args \\
        \$INDEX \\
        $reads \\
        | samtools sort -@ ${task.cpus} ${fasta} -O cram -o ${prefix}.cram -
    """
}