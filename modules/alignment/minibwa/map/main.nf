process MINIBWA_MAP {
    // https://github.com/lh3/minibwa
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/minibwa_samtools_htslib:6f37dc94f6ac9e37'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: aligned
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index, optional: true
    tuple val("${task.process}"), val('minibwa'), eval('minibwa version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_minibwa, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    """
    INDEX=`find -L ./ -name "*.l2b" | sed 's/\\.l2b\$//'`

    minibwa \\
        map \\
        -R "@RG\\tID:${meta.id}\\t${CN}SM:${meta.id}\\tLB:${meta.id}\\tPL:${params.seq_platform}" \\
        ${args} \\
        -t ${task.cpus} \\
        \$INDEX \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -O bam -o ${prefix}.bam -
    """
}