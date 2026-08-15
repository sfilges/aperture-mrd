// CNVkit batch runs access -> autobin -> target/antitarget -> coverage ->
// reference -> fix -> segment -> call in one step. Which of those actually run
// depends on the inputs. Five input scenarios collapse onto four invocations:
//
//   input scenario        tumor  matched normal  reference   flags
//   ----------------------------------------------------------------------------------------
//   tumor/normal          yes    used            -           -f -t -g --output-reference -n <normal>
//   tumor/normal + PoN    yes    DROPPED         pon.cnn     -r pon.cnn
//   tumor-only + PoN      yes    -               pon.cnn     -r pon.cnn
//   tumor-only flat       yes    -               -           -f -t -g --output-reference -n
//   PoN build             -      many, pooled    -           -f -t -g --output-reference -n <n1..nN>
//
// Rows 2 and 3 issue an identical command: when a pooled reference is configured
// the subworkflow discards the matched normal, so a tumor/normal pair is handled
// exactly as tumor-only. That is intentional — a panel built from many
// process-matched normals models assay bias better than a single normal — but it
// does mean a configured --cnvkit_reference silently overrides every matched
// normal in the run, so it is a per-run choice, not a per-sample one.
//
// The combination cannot be pushed down into cnvkit regardless: it exits 1 when
// -n, -t or -g is passed alongside -r, since the reference already fixes the
// bins and bias columns.
//
// Inputs are assumed to be BAM. CNVkit reads CRAM natively but is prohibitively
// slow doing so, and the preprocessing pipeline converts to BAM ahead of any
// caller that needs it.
process CNVKIT_BATCH {
    tag "${meta.id}"
    label 'process_medium'

    container 'community.wave.seqera.io/library/cnvkit_htslib_samtools:86928c121163aca7'

    input:
    tuple val(meta), path(tumor), path(tumor_index), path(normals), path(normal_indices)
    tuple val(meta2), path(fasta), path(fasta_fai)
    tuple val(meta3), path(targets)
    tuple val(meta4), path(access)
    tuple val(meta5), path(reference)
    val method

    output:
    tuple val(meta), path("${prefix}.cnr")            , emit: cnr           , optional: true
    tuple val(meta), path("${prefix}.cns")            , emit: cns           , optional: true
    tuple val(meta), path("${prefix}.call.cns")       , emit: call_cns      , optional: true
    tuple val(meta), path("${prefix}.bintest.cns")    , emit: bintest_cns   , optional: true
    tuple val(meta), path("*.targetcoverage.cnn")     , emit: target_cnn    , optional: true
    tuple val(meta), path("*.antitargetcoverage.cnn") , emit: antitarget_cnn, optional: true
    tuple val(meta), path("${prefix}.reference.cnn")  , emit: reference_cnn , optional: true
    tuple val(meta), path("*.target.bed")             , emit: target_bed    , optional: true
    tuple val(meta), path("*.antitarget.bed")         , emit: antitarget_bed, optional: true
    tuple val(meta), path("*-scatter.pdf")            , emit: scatter       , optional: true
    tuple val(meta), path("*-diagram.pdf")            , emit: diagram       , optional: true
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // `path(normals)` is a single file for a paired normal and a list when
    // building a panel, so normalise before testing or joining.
    def normal_list = normals ? [normals].flatten() : []
    def has_reference = reference ? true : false

    // A pre-built reference already carries the targets, antitargets and bias
    // columns, and cnvkit exits 1 if -t, -n or -g is passed alongside -r (-f is
    // tolerated, but is pointless there). Otherwise build one and keep it: it is
    // the panel in PoN mode and useful provenance elsewhere.
    def reference_args = has_reference
        ? "--reference ${reference}"
        : [
            "--fasta ${fasta}",
            // wgs derives its own targets from the accessible regions
            targets && method != 'wgs' ? "--targets ${targets}" : '',
            access ? "--access ${access}" : '',
            "--output-reference ${prefix}.reference.cnn",
        ].minus('').join(' ')

    // -n/--normal takes nargs='*'. With no filenames it builds a flat reference;
    // with filenames it builds a paired or pooled one. Because the tumor is a
    // positional argument and can be variable --normal must come last on the
    // line or it will swallow the tumor. Not passed at all when reusing -r.
    def normal_args = has_reference ? '' : "--normal ${normal_list.join(' ')}"
    """
    cnvkit.py \\
        batch \\
        ${tumor ?: ''} \\
        --method ${method} \\
        --processes ${task.cpus} \\
        ${reference_args} \\
        ${args} \\
        ${normal_args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Panel mode yields no per-sample segments — only the pooled reference and
    // one coverage pair per normal. Mirroring that keeps stub runs meaningful
    // for the reference-building path.
    def stub_normals = normals ? [normals].flatten() : []
    def per_sample = tumor
        ? """touch ${prefix}.cnr
    touch ${prefix}.cns
    touch ${prefix}.call.cns
    touch ${prefix}.bintest.cns
    touch ${prefix}.targetcoverage.cnn
    touch ${prefix}.antitargetcoverage.cnn"""
        : stub_normals
            .collect { n -> "touch ${n.baseName}.targetcoverage.cnn\n    touch ${n.baseName}.antitargetcoverage.cnn" }
            .join('\n    ')
    // A reused reference is not rewritten
    def stub_reference = reference ? '' : "touch ${prefix}.reference.cnn"
    """
    ${per_sample}
    ${stub_reference}
    """
}
