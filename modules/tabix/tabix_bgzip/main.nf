process TABIX_BGZIPTABIX {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/htslib:1.21--h5efdd21_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"), path("*.gz.tbi"), emit: gz_tbi
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | sed -n 's/^Version: //p'"), emit: versions_tabix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${input.baseName}"
    // Keep the input extension in the output name (e.g. .bed.gz) so tabix can auto-detect
    // the format. Without it tabix falls back to its generic preset, fails to parse every
    // line, and still writes an index (with zero contigs) and exits 0.
    def suffix = input.extension ? ".${input.extension}" : ''

    """
    bgzip --threads ${task.cpus} $args -c ${input} > ${prefix}${suffix}.gz
    tabix $args2 ${prefix}${suffix}.gz

    if [ "\$(tabix -l ${prefix}${suffix}.gz | wc -l)" -eq 0 ]; then
        echo "ERROR: tabix indexed no contigs from ${prefix}${suffix}.gz -- wrong format preset?" >&2
        exit 1
    fi
    """
}
