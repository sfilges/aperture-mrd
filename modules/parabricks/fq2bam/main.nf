process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"

    publishDir "${params.outdir}/${workflow.runName}/preprocessing/parabricks_fq2bam/${meta.id}/", mode: params.publish_dir_mode, pattern: "*.{bam,bai}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.bam")                , emit: bam
    tuple val(meta), path("*.bai")                , emit: bai
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    """
    pbrun fq2bam \\
        --ref $fasta \\
        --in-fq $reads "@RG\\tID:${meta.id}\\t${CN}SM:${meta.id}\\tLB:${meta.id}\\tPL:${params.seq_platform}" \\
        --out-bam ${prefix}.bam \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}