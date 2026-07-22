process FASTP {
    tag "$meta.id with $task.cpus cores"
    label 'process_medium'

    container 'biocontainers/fastp:1.3.6--h43da1c4_0'

    input:
    tuple val(meta), path(reads)
    val   use_merged

    output:
    tuple val(meta), path('*_fastp.fastq.gz')   , emit: reads
    tuple val(meta), path('*.json')             , emit: json
    tuple val(meta), path('*.html')             , emit: html
    tuple val(meta), path('*.log')              , emit: log
    tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | sed -e "s/fastp //g"'), emit: versions_fastp, topic: versions
    tuple val(meta), path('*_merged.fastq.gz'), optional:true, emit: reads_merged

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_fastq = use_merged ? "-m --merged_out ${prefix}_merged.fastq.gz" : ''
    """
    [ ! -f ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    fastp \\
        --in1 ${prefix}_1.fastq.gz \\
        --in2 ${prefix}_2.fastq.gz \\
        --out1 ${prefix}_1_fastp.fastq.gz \\
        --out2 ${prefix}_2_fastp.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        ${merge_fastq} \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        ${args} \\
        2> >(tee ${prefix}.fastp.log >&2)
    """
}