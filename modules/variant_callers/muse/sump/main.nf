process MUSE_SUMP {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/muse:6637291dcbb0bdb8'

    input:
    tuple val(meta), path(call_txt)
    path dbsnp                      // dbSNP vcf file that should be bgzip compressed,tabix indexed and based on the same reference genome used in 'MuSE call'
    path dbsnp_tbi
    val error_model

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('MuSE'), eval('MuSE --version | sed -e "s/MuSE, version //g" | sed -e "s/MuSE v//g"'), emit: versions_muse, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model = error_model == "wgs" ? "-G" : "-E"
    """
    MuSE sump \\
        -I $call_txt \\
        $model \\
        -D $dbsnp \\
        -O "${prefix}.muse.vcf" \\
        -n $task.cpus \\
        $args        
    """
}