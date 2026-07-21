process GATK4_MARKDUPLICATES {
    tag "$meta.id with $task.cpus cores"
    label 'process_medium'

    container 'biocontainers/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:7cc3d06cbf42e28c5e2ebfc7c858654c7340a9d5-0'

    publishDir "${params.outdir}/${workflow.runName}/preprocessing/markduplicates/${meta.id}/", mode: params.publish_dir_mode, pattern: "*.{cram,crai}"
    publishDir "${params.outdir}/${workflow.runName}/reports/", mode: params.publish_dir_mode, pattern: "*.metrics"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.cram"),    emit: cram,  optional: true
    tuple val(meta), path("*.crai"),    emit: crai,  optional: true
    tuple val(meta), path("*.metrics"), emit: metrics
    tuple val("${task.process}"), val('gatk4'), eval("echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//'"), emit: versions_gatk4, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//'"), emit: versions_samtools, topic: versions

    script:
    def input_list = bam.collect{"--INPUT $it"}.join(' ')
    def prefix = "${meta.id}"
    def outfile = "${prefix}.sorted.md.bam"

    // Using samtools and not Markduplicates to compress to CRAM speeds up computation:
    // https://medium.com/@acarroll.dna/looking-at-trade-offs-in-compression-levels-for-genomics-tools-eec2834e8b94
    """
    gatk --java-options "-Xms4g -Xmx8g" \\
        MarkDuplicates \\
        $input_list \\
        --OUTPUT ${outfile} \\
        --METRICS_FILE ${prefix}.md.metrics

    samtools view --cram --with-header --threads ${task.cpu} --reference ${fasta} --output ${prefix}.sorted.md.cram ${outfile}
    rm ${outfile}
    samtools index ${prefix}.sorted.md.cram 
    """
}