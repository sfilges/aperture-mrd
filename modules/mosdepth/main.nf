process MOSDEPTH {
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/mosdepth:0.3.11--h0ec343a_1'

    publishDir "${params.outdir}/${workflow.runName}/reports/mosdepth/${meta.id}", mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(cram), path(crai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    script:
    def prefix = "${meta.id}"

    """
    mosdepth \\
        --threads $task.cpus \\
        -n --fast-mode --by 500 \\
        --fasta ${fasta} \\
        $prefix \\
        $cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}