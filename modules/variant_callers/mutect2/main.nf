process GATK4_MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'
    
    publishDir "${params.outdir}/${workflow.runName}/variant_calling/", mode: params.publish_dir_mode, pattern: "*{vcf.gz,vcf.gz.tbi,stats}", saveAs: { meta.num_intervals > 1 ? null : "mutect2/${meta.id}/${it}" }

    // Input format:
    //[meta, normal_cram, normal_crai, tumor_cram, tumor_crai]

    input:
    tuple val(meta), path(normal_cram), path(normal_crai), path(tumor_cram), path(tumor_crai), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(dict)
    path(germline_resource)
    path(germline_resource_tbi)
    path(panel_of_normals)
    path(panel_of_normals_tbi)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    tuple val("${task.process}"), val('gatk4'), eval("echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//'"), emit: versions_gatk4, topic: versions

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def gr_command       = germline_resource ? "--germline-resource $germline_resource" : ""
    def pon_command      = panel_of_normals ? "--panel-of-normals $panel_of_normals" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        --input $tumor_cram \\
        --input $normal_cram \\
        --output ${prefix}.vcf.gz \\
        --normal-sample ${meta.normal_id} \\
        --reference $fasta \\
        $interval_command \\
        ${gr_command} \\
        ${pon_command} \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        --dont-use-soft-clipped-bases \\
        --tmp-dir . \\
        $args
    """
}