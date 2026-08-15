process FULCRUM_FGUMI_ALIGN {
    // https://github.com/fulcrumgenomics/fgumi
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa-mem3_fgumi:50f5ff04bb8a3e5d'

    input:
    tuple val(meta), path(bam)
    path(fasta)

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''
    def args3 = task.ext.args ?: ''
    def args4 = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // --max-memory specifies memory per thread. Total memory = `max_memory` × threads
    """
    fgumi fastq --input -t ${task.cpus} ${args} ${bam} \\
        | bwa-mem3 -K 100000000 -Y --bam=0 -p -t ${task.cpus} ${args2} ${fasta} - \\
        | fgumi zipper --bwa-chunk-size 100000000 --threads ${task.cpus} ${args3} --unmapped ${bam} -r ${fasta} \\
        | fgumi sort --max-memory 1G --threads ${task.cpus} ${args4} --output ${prefix}.sorted.bam --order template-coordinate
    """
}