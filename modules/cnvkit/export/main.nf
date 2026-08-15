// Export CNVkit segments as DNAcopy-format .seg for PureCN.
//
// PureCN consumes CNVkit output via --seg-file (this output) plus --tumor
// (the .cnr), which is why it needs no interval file of its own. The PureCN
// docs pair this with --enumerate-chroms; that is set via ext.args so it stays
// visible in the config rather than buried here.
process CNVKIT_EXPORT_SEG {
    tag "${meta.id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/cnvkit_htslib_samtools:86928c121163aca7'

    input:
    tuple val(meta), path(cns)

    output:
    tuple val(meta), path("*.seg"), emit: seg
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cnvkit.py \\
        export seg \\
        ${cns} \\
        ${args} \\
        --output ${prefix}.seg
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.seg
    """
}
