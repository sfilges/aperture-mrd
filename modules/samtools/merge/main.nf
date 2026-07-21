process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(input_files, stageAs: "?/*")
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    tuple val(meta), path("*.csi")         , optional:true, emit: csi
    tuple val(meta), path("*.crai")        , optional:true, emit: crai
    tuple val("${task.process}"), val('samtools'), eval("echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//'"), emit: versions_samtools, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${reference} \\
        ${prefix}.${file_type} \\
        $input_files
    """
}