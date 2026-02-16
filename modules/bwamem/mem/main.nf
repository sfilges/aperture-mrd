process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/bwa_htslib_samtools:fea458eb782fa887'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = { params.split_fastq > 1 ? "${meta.id}".concat('.').concat(reads.get(0).name.tokenize('.')[0]) : "${meta.id}.sorted" }()
    def extra_args = meta.status == 1 ? "-K 100000000 -Y -B 3" : "-K 100000000 -Y"
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa \\
        mem \\
        -R "@RG\\tID:${meta.id}\\t${CN}SM:${meta.id}\\tLB:${meta.id}\\tPL:${params.seq_platform}" \\
        -t $task.cpus \\
        $extra_args \\
        \$INDEX \\
        $reads \\
        | samtools sort -@ $task.cpus -O bam -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
