process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'

    // TODO Use a container with gatk and samtools installed, to index the output CRAM file directly
    container 'biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' // Consider quay.io/nf-core/gatk:4.6.1.10

    publishDir "${params.outdir}/${workflow.runName}/preprocessing/recal/${meta.id}/", mode: params.publish_dir_mode, pattern: "*.{cram,crai}"

    input:
    tuple val(meta), path(cram), path(bqsr_table)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path  dict

    output:
    tuple val(meta), path("*.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//'"), emit: versions_gatk4, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.recal"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input $cram \\
        --output ${prefix}.${cram.getExtension()} \\
        --reference $fasta \\
        --bqsr-recal-file $bqsr_table \\
        --tmp-dir . \\
        $args
    """
}