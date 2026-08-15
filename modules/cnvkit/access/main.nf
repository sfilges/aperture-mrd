process CNVKIT_ACCESS {
    tag "${meta.id}"
    label 'process_low'

    container 'community.wave.seqera.io/library/cnvkit_htslib_samtools:86928c121163aca7'

    input:
    tuple val(meta), path(fasta)
    path exclude_beds

    output:
    tuple val(meta), path("*.access.bed"), emit: bed
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    // Other known unmappable, variable, or poorly sequenced regions can be
    // excluded with -x/--exclude, which is repeatable, so emit one flag per BED.
    // Wrapping in a list first stops a single unstaged file from being iterated
    // as path components rather than treated as one file.
    def exclude = exclude_beds
        ? [exclude_beds].flatten().collect { bed -> "--exclude ${bed}" }.join(' ')
        : ''
    """
    cnvkit.py \\
        access \\
        ${fasta} \\
        ${exclude} \\
        ${args} \\
        --output ${prefix}.access.bed
    """

    stub:
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    touch ${prefix}.access.bed
    """
}
