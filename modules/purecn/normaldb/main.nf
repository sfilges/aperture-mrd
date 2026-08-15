process PURECN_NORMALDB {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/bioconductor-dnacopy_bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_pruned:ca4b5595ad5ac8ff'

    input:
    tuple val(meta), path(coverage_files), path(normal_vcf), path(normal_vcf_tbi)
    val   genome
    val   assay

    output:
    tuple val(meta), path("normalDB*.rds")               , emit: rds
    tuple val(meta), path("interval_weights*.png")       , emit: png
    tuple val(meta), path("mapping_bias*.rds")           , emit: bias_rds,    optional: true
    tuple val(meta), path("mapping_bias_hq_sites*.bed")  , emit: bias_bed,    optional: true
    tuple val(meta), path("low_coverage_targets*.bed")   , emit: low_cov_bed, optional: true
    tuple val("${task.process}"), val('purecn'), eval("Rscript -e 'cat(as.character(packageVersion(\"PureCN\")))'"), emit: versions_purecn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args     ?: ''
    def normal_panel    = normal_vcf        ? "--normal-panel ${normal_vcf}" : ""
    """
    echo $coverage_files | tr ' ' '\\n' > coverages.list
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/NormalDB.R --out-dir ./ \\
        --coverage-files coverages.list \\
        --genome ${genome} \\
        --assay ${assay} \\
        ${normal_panel} \\
        $args
    """

    stub:
    // NormalDB.R names its outputs from --assay and --genome, not from the meta
    // id, and the mapping bias files appear only when a --normal-panel VCF was
    // given. That flag is built from the normal_vcf input rather than ext.args,
    // so the stub keys off the same input the script does.
    def mapping_bias = normal_vcf ? "touch mapping_bias_${assay}_${genome}.rds" : ''
    def mapping_bias_hq_sites = normal_vcf ? "touch mapping_bias_hq_sites_${assay}_${genome}.bed" : ''
    """
    touch normalDB_${assay}_${genome}.rds
    touch interval_weights_${assay}_${genome}.png
    touch low_coverage_targets_${assay}_${genome}.bed
    ${mapping_bias}
    ${mapping_bias_hq_sites}
    """
}