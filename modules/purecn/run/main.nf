// Estimate tumor purity and ploidy from third-party (CNVkit) copy number data.
//
// Fed the CNVkit .cnr as --tumor and the exported .seg as --seg-file, PureCN
// needs no interval file and no coverage files of its own: the .cnr already
// carries the intervals. --intervals is therefore optional here and only used
// when driving PureCN from its own IntervalFile.R output instead.
//
// The VCF must contain germline as well as somatic sites, which is why Mutect2
// is run with --genotype-germline-sites for this path.
process PURECN_RUN {
    tag "${meta.id}"
    label 'process_medium'

    container 'community.wave.seqera.io/library/bioconductor-dnacopy_bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_pruned:ca4b5595ad5ac8ff'

    input:
    tuple val(meta), path(tumor), path(seg_file), path(vcf), path(vcf_tbi)
    tuple val(meta2), path(intervals)
    path normal_db
    path mapping_bias
    path snp_blacklist
    val genome

    output:
    tuple val(meta), path("${prefix}.csv")                 , emit: csv
    tuple val(meta), path("${prefix}.rds")                 , emit: rds
    tuple val(meta), path("${prefix}.pdf")                 , emit: pdf
    tuple val(meta), path("${prefix}_dnacopy.seg")         , emit: seg
    tuple val(meta), path("${prefix}_local_optima.pdf")    , emit: local_optima_pdf        , optional: true
    tuple val(meta), path("${prefix}_genes.csv")           , emit: genes_csv               , optional: true
    tuple val(meta), path("${prefix}_amplification_pvalues.csv"), emit: amplification_pvalues_csv, optional: true
    tuple val(meta), path("${prefix}.vcf.gz")              , emit: vcf_gz                  , optional: true
    tuple val(meta), path("${prefix}_variants.csv")        , emit: variants_csv            , optional: true
    tuple val(meta), path("${prefix}_loh.csv")             , emit: loh_csv                 , optional: true
    tuple val(meta), path("${prefix}_chromosomes.pdf")     , emit: chr_pdf                 , optional: true
    tuple val(meta), path("${prefix}_segmentation.pdf")    , emit: segmentation_pdf        , optional: true
    tuple val(meta), path("${prefix}_multisample.seg")     , emit: multisample_seg         , optional: true
    tuple val(meta), path("${prefix}.log")                 , emit: log                     , optional: true
    tuple val("${task.process}"), val('purecn'), eval("Rscript -e 'cat(as.character(packageVersion(\"PureCN\")))'"), emit: versions_purecn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def seg_opt = seg_file ? "--seg-file ${seg_file}" : ''
    def intervals_opt = intervals ? "--intervals ${intervals}" : ''
    def vcf_opt = vcf ? "--vcf ${vcf}" : ''
    def normaldb_opt = normal_db ? "--normaldb ${normal_db}" : ''
    def mapping_bias_opt = mapping_bias ? "--mapping-bias-file ${mapping_bias}" : ''
    // Only worth setting when neither a matched normal nor a large pool of
    // normals is available; the PureCN docs recommend the UCSC simple repeats
    // track here. Note this filters SNV positions, unlike the CNVkit access
    // exclusions which remove regions from coverage binning.
    def blacklist_opt = snp_blacklist ? "--snp-blacklist ${snp_blacklist}" : ''
    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/PureCN.R \\
        --out ./ \\
        --sampleid ${prefix} \\
        --tumor ${tumor} \\
        --genome ${genome} \\
        --parallel \\
        --cores ${task.cpus} \\
        ${seg_opt} \\
        ${intervals_opt} \\
        ${normaldb_opt} \\
        ${mapping_bias_opt} \\
        ${blacklist_opt} \\
        ${vcf_opt} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // The curation CSV is parsed downstream for purity and ploidy, so the stub
    // emits a real row rather than an empty file — otherwise a stub run silently
    // produces an empty channel instead of exercising the path.
    """
    cat <<-END_CURATION > ${prefix}.csv
    Sampleid,Purity,Ploidy,Sex,Contamination,Flagged,Failed,Curated,Comment
    ${prefix},0.65,2.1,?,0,FALSE,FALSE,FALSE,
    END_CURATION

    touch ${prefix}.rds
    touch ${prefix}.pdf
    touch ${prefix}_dnacopy.seg
    touch ${prefix}_local_optima.pdf
    """
}
