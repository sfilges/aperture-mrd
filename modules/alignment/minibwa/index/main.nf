process MINIBWA_INDEX {
    tag "$fasta"
    label 'process_medium'
    // NOTE minibwa builds an FM-index with libsais; peak memory scales with the reference size.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    container 'community.wave.seqera.io/library/minibwa:0.7--8e8120c7ca8465fb'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("minibwa"), emit: index
    tuple val("${task.process}"), val('minibwa'), eval('minibwa version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_minibwa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir minibwa
    minibwa \\
        index \\
        $args \\
        -t ${task.cpus} \\
        $fasta \\
        minibwa/${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    // MINIBWA_MAP globs for *.l2b to find the index prefix.
    """
    mkdir minibwa
    touch minibwa/${prefix}.l2b
    touch minibwa/${prefix}.bwt
    touch minibwa/${prefix}.ann
    """
}