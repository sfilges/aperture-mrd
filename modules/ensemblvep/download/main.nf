process ENSEMBLVEP_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/ensembl-vep:116.0--pl5321h2a3209d_0'

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)
    val preflight_check

    output:
    tuple val(meta), path(prefix), emit: cache
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help 2>&1 | sed -n 's/^ *ensembl-vep *: *//p'"), emit: versions_ensemblvep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'vep_cache'
    // The cache tarball is ~20 GB and vep_install only discovers a bad
    // species/assembly/version combination once it is already downloading.
    // Checking the release manifest first turns an hours-long failure into a
    // seconds-long one. HTTP::Tiny is core Perl, so no extra dependency.
    def filename = "${species}_vep_${cache_version}_${assembly}.tar.gz"
    def checksums_url = "https://ftp.ensembl.org/pub/release-${cache_version}/variation/indexed_vep_cache/CHECKSUMS"
    """
    if [ "${preflight_check}" = "true" ]; then
        perl -MHTTP::Tiny -e '
            my \$r = HTTP::Tiny->new(timeout => 30)->get("${checksums_url}");
            \$r->{success} or die "Failed to fetch CHECKSUMS (HTTP \$r->{status})\\n";
            \$r->{content} =~ /\\Q${filename}\\E/ or die "${filename} not found in CHECKSUMS\\n";
            print "Pre-flight OK: ${filename} found in CHECKSUMS\\n";
        '
    fi

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
    mkdir -p ${prefix}/${species}/${cache_version}_${assembly}
    """
}
