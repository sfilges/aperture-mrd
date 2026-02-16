process MUSE_CALL {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/muse:6637291dcbb0bdb8'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_index)

    output:
    tuple val(meta), path("*.MuSE.txt"), emit: txt
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MuSE \\
        call \\
        $args \\
        -f $fasta \\
        -O ${prefix}  \\
        -n $task.cpus \\
        $tumor_bam    \\
        $normal_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" | sed -e "s/MuSE v//g")
    END_VERSIONS
    """
}