process FULCRUM_FGUMI_FILTER {
    tag {meta.id}
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa-mem3_fgumi:50f5ff04bb8a3e5d'

    input:
    tuple val(meta), path(consensus_bam)
    path(index)
    path(fasta)

    output:
    tuple val(meta), path("*filtered.bam"), emit: bam, optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fgumi fastq --input ${consensus_bam} ${args} \
        | bwa-mem3 -t ${task.cpus} -p -K 150000000 -Y ${index} - \
        | fgumi zipper --unmapped ${consensus_bam} --reference ${fasta} \
        | fgumi filter--ref ${fasta} --min-reads 3 \
        | fgumi sort --output ${prefix}.filtered.bam --order coordinate --threads 4
    """
}