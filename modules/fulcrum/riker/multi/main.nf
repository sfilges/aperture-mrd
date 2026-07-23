process RIKER_MULTI {
    // https://github.com/fulcrumgenomics/riker
    // As a rough guide returns flatten past ~6 threads for BAM and ~8 for CRAM
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/riker:0.4.1--hec9b1f2_0'

    input:
    tuple val(meta), path(input) // bam or cram
    tuple val(meta2), path(fasta) // Reference FASTA file (must be indexed with .fai). Required for CRAM and some tools.
    tuple val(meta3), path(fai)
    tuple val(meta4), path(baits), path(targets)

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
    def ref = fasta ?: ''
    def hybcap_opts = (baits && targets) ? "--hybcap::baits ${baits} --hybcap::targets ${targets}" : ''
    """
    riker multi \\
        -i ${input} \\
        -r ${ref} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${hybcap_opts}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment-metrics.txt
    touch ${prefix}.base-distribution-by-cycle.txt
    touch ${prefix}.mean-quality-by-cycle.txt
    touch ${prefix}.quality-score-distribution.txt
    touch ${prefix}.error-mismatch.txt
    touch ${prefix}.error-overlap.txt
    touch ${prefix}.error-indel.txt
    touch ${prefix}.gcbias-detail.txt
    touch ${prefix}.gcbias-summary.txt
    touch ${prefix}.hybcap-metrics.txt
    touch ${prefix}.hybcap-per-target.txt
    touch ${prefix}.hybcap-per-base.txt
    touch ${prefix}.isize-metrics.txt
    touch ${prefix}.isize-histogram.txt
    touch ${prefix}.wgs-metrics.txt
    touch ${prefix}.wgs-coverage.txt
    touch ${prefix}.base-distribution-by-cycle.pdf
    touch ${prefix}.gcbias-chart.pdf
    touch ${prefix}.isize-histogram.pdf
    touch ${prefix}.mean-quality-by-cycle.pdf
    touch ${prefix}.quality-score-distribution.pdf
    touch ${prefix}.wgs-coverage.pdf
    """
}