process GATK4_COLLECTWGSMETRICS {
    tag "$meta.sample"

    label 'gatk4'

    input:
    tuple val(meta), path(reads)
    path fasta
    path index

    output:
    path "metrics/*", emit: metrics
    path "versions.txt", emit: versions

    when:
    params.collectwgsmetrics

    script:
    def cmd = "gatk CollectWgsMetrics \\
        -I ${reads} \\
        -R ${fasta} \\
        --BAIT_INTERVALS ${index} \\
        --TARGET_INTERVALS ${index} \\
        -O metrics/${meta.sample}.wgs_metrics.txt"

    return cmd
}