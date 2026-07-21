process TABIX_BGZIPTABIX {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/htslib:1.21--h5efdd21_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"), path("*.gz.tbi"), emit: gz_tbi
    tuple val("${task.process}"), val('tabix'), eval("echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//'"), emit: versions_tabix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${input.baseName}"

    """
    bgzip --threads ${task.cpus} $args -c ${input} > ${prefix}.gz
    tabix $args2 ${prefix}.gz
    """
}
