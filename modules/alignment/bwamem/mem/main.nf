process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa_htslib_samtools:fea458eb782fa887'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val("${task.process}"), val('bwa'), eval("echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact.*\$//'"), emit: versions_bwa, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = { params.split_fastq > 1 ? "${meta.id}".concat('.').concat(reads.get(0).name.tokenize('.')[0]) : "${meta.id}.sorted" }()
    def extra_args = meta.status == 1 ? "-K 100000000 -Y -B 3" : "-K 100000000 -Y"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa \\
        mem \\
        -R "${meta.read_group}" \\
        -t $task.cpus \\
        $extra_args \\
        \$INDEX \\
        $reads \\
        | samtools sort -@ $task.cpus -O bam -o ${prefix}.bam -
    """
}
