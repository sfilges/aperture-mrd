process NGSCHECKMATE {
    label 'process_low'

    container 'biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1'

    input:
    tuple val(meta) , path(files)
    tuple val(meta2), path(snp_bed)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*_corr_matrix.txt"), emit: corr_matrix
    tuple val(meta), path("*_matched.txt")    , emit: matched
    tuple val(meta), path("*_all.txt")        , emit: all
    tuple val(meta), path("*.pdf")            , emit: pdf, optional: true
    tuple val(meta), path("*.vcf")            , emit: vcf, optional: true
    tuple val("${task.process}"), val('ngscheckmate'), eval('ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g"'), emit: versions_ngscheckmate, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    def unzip = files.any { it.toString().endsWith(".vcf.gz") }
    """
    if $unzip
    then
        for VCFGZ in *.vcf.gz; do
            gunzip -cdf \$VCFGZ > \$( basename \$VCFGZ .gz );
        done
    fi

    NCM_REF="./"${fasta} ncm.py -d . -bed ${snp_bed} -O . -N ${prefix} $args

    if $unzip
    then
        rm -f *.vcf  # clean up decompressed vcfs
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    touch ${prefix}_output_corr_matrix.txt
    touch ${prefix}_matched.txt
    touch ${prefix}_all.txt
    touch ${prefix}.pdf
    """

}