process ENSEMBLVEP_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/ensembl-vep:116.0--pl5321h2a3209d_0'

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path(prefix), emit: cache
    tuple val("${task.process}"), val('ensemblvep'), eval("echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//'"), emit: versions_ensemblvep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    vep_install \\
        --CACHEDIR $prefix \\
        --SPECIES $species \\
        --ASSEMBLY $assembly \\
        --CACHE_VERSION $cache_version \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    mkdir $prefix
    """
}