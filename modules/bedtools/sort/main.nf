process BEDTOOLS_SORT {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    tuple val(meta), path(intervals)
    path genome_file

    output:
    tuple val(meta), path("*.sorted.bed"), emit: sorted
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def genome_cmd = genome_file     ?  "-g $genome_file" : ""
    """
    bedtools sort \\
        -i $intervals \\
        $genome_cmd \\
        $args \\
        > "${prefix}.sorted.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}