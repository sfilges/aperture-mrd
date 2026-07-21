process FASTP {
    tag "$meta.id with $task.cpus cores"
    label 'process_medium'

    container 'biocontainers/fastp:0.23.4--hadf994f_2'

    publishDir "${params.outdir}/${workflow.runName}/reports/fastp/${meta.id}", mode: params.publish_dir_mode, pattern: "*.{html,json,log}"
    publishDir "${params.outdir}/${workflow.runName}/preprocessing/fastp/${meta.id}", mode: params.publish_dir_mode, pattern: "*.fastq.gz", enabled: params.save_fastqs

    input:
    tuple val(meta), path(reads)
    val   use_merged
    val   split_fastq
    val   trim_nextseq
    val   length_required
    val   save_fastq

    output:
    tuple val(meta), path('*_fastp.fastq.gz')   , emit: reads
    tuple val(meta), path('*.json')             , emit: json
    tuple val(meta), path('*.html')             , emit: html
    tuple val(meta), path('*.log')              , emit: log
    tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | sed -e "s/fastp //g"'), emit: versions_fastp, topic: versions
    tuple val(meta), path('*_merged.fastq.gz'), optional:true, emit: reads_merged

    script:
    def prefix = "${meta.id}"
    def merge_fastq = use_merged ? "-m --merged_out ${prefix}_merged.fastq.gz" : ''
    // if number of lines < split_fastq, fastp will create at least one chunk per cpu allocated
    def chunk_fastq = split_fastq > 0 ? "--split_by_lines ${params.split_fastq * 4}" : ''
    def trim_poly_g = trim_nextseq ? '--trim_poly_g' : ''
    def min_length = length_required > 0 ? "--length_required ${params.length_required}": ''
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
        $min_length \\
        $merge_fastq \\
        ${chunk_fastq} \\
        $trim_poly_g \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        $args \\
        2> >(tee ${prefix}.fastp.log >&2)
    """
}