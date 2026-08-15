// Convert segmented log2 ratios to absolute integer copy number.
//
// Given a purity estimate (from PureCN) the 'clonal' method rescales by purity
// and ploidy before rounding, which is materially more accurate than the fixed
// log2 cutoffs. Without one it falls back to 'threshold', so the process still
// produces usable calls when PureCN is disabled or fails to converge.
process CNVKIT_CALL {
    tag "${meta.id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/cnvkit_htslib_samtools:86928c121163aca7'

    input:
    tuple val(meta), path(cns), val(purity), val(ploidy)

    output:
    tuple val(meta), path("*.called.cns"), emit: cns
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def method_args = purity ? "--method clonal --purity ${purity}" : "--method threshold"
    // Only meaningful alongside --purity; cnvkit defaults to 2 otherwise.
    def ploidy_args = purity && ploidy ? "--ploidy ${ploidy}" : ''
    """
    cnvkit.py \\
        call \\
        ${cns} \\
        ${method_args} \\
        ${ploidy_args} \\
        ${args} \\
        --output ${prefix}.called.cns
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.called.cns
    """
}
