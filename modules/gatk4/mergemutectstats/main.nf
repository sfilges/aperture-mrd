process GATK4_MERGEMUTECTSTATS {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*.merged.stats"), emit: stats
    tuple val("${task.process}"), val('gatk4'), eval("echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//'"), emit: versions_gatk4, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = stats.collect { "--stats $it" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MergeMutectStats] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MergeMutectStats \\
        $input_list \\
        --output ${prefix}.merged.stats \\
        --tmp-dir . \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.stats
    """
}
