process GATK4_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'

    // Requires GATK and R to be installed in the container
    container 'quay.io/nf-core/gatk:4.6.1.0'

    publishDir "${params.outdir}/${workflow.runName}/reports/gatk/${meta.id}", mode: params.publish_dir_mode, pattern: "*.pdf"

    input:
    tuple val(meta), path(cram)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    path "*.alignment_summary_metrics", emit: alignment_summary_metrics, optional: true
    path "*.insert_size_metrics", emit: insert_size_metrics, optional: true
    path "*.quality_by_cycle_metrics", emit: quality_by_cycle_metrics, optional: true
    path "*.base_distribution_by_cycle_metrics", emit: base_distribution_by_cycle_metrics, optional: true
    path "*.gc_bias.summary_metrics", emit: gc_bias_metrics, optional: true
    path "*.quality_yield_metrics", emit: quality_yield_metrics, optional: true
    path "versions.yml", emit: versions

    script:
    """
    gatk CollectMultipleMetrics \\
        -I ${cram} \\
        -R ${fasta} \\
        -O ${meta.id} \\
        --PROGRAM CollectAlignmentSummaryMetrics \\
        --PROGRAM CollectInsertSizeMetrics \\
        --PROGRAM MeanQualityByCycle \\
        --PROGRAM CollectBaseDistributionByCycle \\
        --PROGRAM CollectGcBiasMetrics \\
        --PROGRAM CollectSequencingArtifactMetrics \\
        --PROGRAM CollectQualityYieldMetrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}