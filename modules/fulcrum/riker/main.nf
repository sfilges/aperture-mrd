process RIKER_MULTI {
    // https://github.com/fulcrumgenomics/riker
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/riker:0.4.1--hec9b1f2_0'

    input:
    tuple val(meta), path(input) // bam or cram
    tuple val(meta2), path(fasta) // Reference FASTA file (must be indexed with .fai). Required for CRAM and some tools.
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.alignment-metrics.txt"),          emit: alignment_metrics,         optional: true
    tuple val(meta), path("*.base-distribution-by-cycle.txt"), emit: base_dist,                 optional: true
    tuple val(meta), path("*.error-indel.txt"),                emit: error_indel,               optional: true
    tuple val(meta), path("*.error-mismatch.txt"),             emit: error_mismatch,            optional: true
    tuple val(meta), path("*.error-overlap.txt"),              emit: error_overlap,             optional: true
    tuple val(meta), path("*.gcbias-detail.txt"),              emit: gcbias_detail,             optional: true
    tuple val(meta), path("*.gcbias-summary.txt"),             emit: gcbias_summary,            optional: true
    tuple val(meta), path("*.hybcap-metrics.txt"),             emit: hybcap_metrics,            optional: true
    tuple val(meta), path("*.hybcap-per-base.txt*"),           emit: hybcap_per_base,           optional: true
    tuple val(meta), path("*.hybcap-per-target.txt"),          emit: hybcap_per_target,         optional: true
    tuple val(meta), path("*.isize-histogram.txt"),            emit: isize_histogram,           optional: true
    tuple val(meta), path("*.isize-metrics.txt"),              emit: isize_metrics,             optional: true
    tuple val(meta), path("*.mean-quality-by-cycle.txt"),      emit: mean_qual,                 optional: true
    tuple val(meta), path("*.pdf"),                            emit: pdf,                       optional: true
    tuple val(meta), path("*.quality-score-distribution.txt"), emit: qual_dist,                 optional: true
    tuple val(meta), path("*.rna-biotype.txt"),                emit: rna_biotype,               optional: true
    tuple val(meta), path("*.rna-insert-size-histogram.txt"),  emit: rna_insert_size_histogram, optional: true
    tuple val(meta), path("*.rna-insert-size.txt"),            emit: rna_insert_size,           optional: true
    tuple val(meta), path("*.rna-metrics.txt"),                emit: rna_metrics,               optional: true
    tuple val(meta), path("*.wgs-coverage.txt"),               emit: wgs_coverage,              optional: true
    tuple val(meta), path("*.wgs-metrics.txt"),                emit: wgs_metrics,               optional: true
    tuple val("${task.process}"), val('riker'), eval("riker --version 2>&1 | sed 's/riker //'") , topic: versions, emit: versions_riker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    riker multi \\
        -i ${input} \\
        -r ${fasta} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        --tools alignment isize basic wgs gcbias alignment error \\
        ${args}
    """
}