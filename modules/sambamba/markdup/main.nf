process SAMBAMBA_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/sambamba:1.0.1--h6f6fda4_0'

    input:
    tuple val(meta), path(bam) // Aligned and sorted BAM in

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.sorted.md.bam"

    """
    sambamba \\
        markdup \\
        $args \\
        -t $task.cpus \\
        --tmpdir ./ \\
        $bam \\
        ${prefix}.bam

    sambamba index \\
        --nthreads $task.cpus \\
        --tmpdir ./ \\
        ${prefix} ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}