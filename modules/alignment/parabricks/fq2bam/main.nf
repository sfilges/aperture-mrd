process PARABRICKS_FQ2BAM {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to run on a cluster because we need to copy the fasta reference, see https://github.com/nf-core/modules/issues/9230
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    tuple val(meta4), path(intervals)
    tuple val(meta5), path(known_sites)
    val output_fmt

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.table"), emit: bqsr_table, optional: true
    tuple val(meta), path("*_qc_metrics"), emit: qc_metrics, optional: true
    tuple val(meta), path("*.duplicate-metrics.txt"), emit: duplicate_metrics, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def in_fq_command = meta.single_end ? "--in-se-fq ${reads}" : "--in-fq ${reads}"
    def extension = "${output_fmt}"

    def known_sites_command = known_sites ? (known_sites instanceof List ? known_sites.collect { knownSite -> "--knownSites ${knownSite}" }.join(' ') : "--knownSites ${known_sites}") : ""
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def intervals_command = intervals ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' ') : "--interval-file ${intervals}") : ""

    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp ${fasta} \$INDEX

    pbrun \\
        fq2bam \\
        --ref \$INDEX \\
        ${in_fq_command} \\
        --out-bam ${prefix}.${extension} \\
        ${known_sites_command} \\
        ${known_sites_output_cmd} \\
        ${intervals_command} \\
        ${num_gpus} \\
        --bwa-cpu-thread-pool ${task.cpus} \\
        --monitor-usage \\
        ${args}
    """
}