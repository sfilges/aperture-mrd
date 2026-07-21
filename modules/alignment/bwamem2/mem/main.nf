process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.bam")  , emit: bam , optional:true
    tuple val(meta), path("*.cram") , emit: cram, optional:true
    tuple val(meta), path("*.crai") , emit: crai, optional:true
    tuple val(meta), path("*.csi")  , emit: csi , optional:true
    tuple val("${task.process}"), val('bwamem2'), eval("echo \$(bwa-mem2 version 2>&1) | sed 's/.* //'"), emit: versions_bwamem2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//'"), emit: versions_samtools, topic: versions

    script:
    def prefix = { params.split_fastq > 1 ? "${meta.id}".concat('.').concat(reads.get(0).name.tokenize('.')[0]) : "${meta.id}.sorted" }()
    def args = meta.status == 1 ? "-K 100000000 -Y -B 3" : "-K 100000000 -Y"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem2 \\
        mem \\
        -R "${meta.read_group}" \\
        -t $task.cpus \\
        $args \\
        \$INDEX \\
        $reads \\
        | samtools sort -@ $task.cpus -O bam -o ${prefix}.bam -
    """
}