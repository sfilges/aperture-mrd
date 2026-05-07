process TMB_PYTMB {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/tmb:1.5.0--0724a7e1e50a32cd'

    input:
    tuple val(meta), path(vcf), path(bed), val(eff_genome_size), path(var_config), path(db_config)

    output:
    tuple val(meta), path("*.log")          , emit: tmb_log
    tuple val(meta), path("*_export.vcf.gz"), optional:true, emit: export_vcf
    tuple val(meta), path("*_debug.vcf.gz") , optional:true, emit: debug_vcf
    tuple val("${task.process}"), val('tmb'), eval(" pyTMB.py --version | awk '{print \$2}' | tr -d '()' "), topic: versions, emit: versions_tmb

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (eff_genome_size && bed) {
        log.warn("Both 'eff_genome_size' and 'bed' were provided. 'eff_genome_size' will be used.")
    } else if (!eff_genome_size && !bed) {
        log.error("One of 'eff_genome_size' or 'bed' must be provided.")
    }

    def genome_size = eff_genome_size ? "--effGenomeSize ${eff_genome_size}" : "--bed ${bed}"

    """
    pyTMB.py \\
        $args \\
        -i $vcf \\
        $genome_size \\
        --dbConfig $db_config \\
        --varConfig $var_config > ${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.log
    echo "" | gzip > ${prefix}_export.vcf.gz
    echo "" | gzip > ${prefix}_debug.vcf.gz
    """
}