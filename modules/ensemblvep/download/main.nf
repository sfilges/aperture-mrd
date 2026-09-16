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
    prefix = task.ext.prefix ?: 'vep_cache'
    // The cache is fetched over HTTPS rather than with vep_install. vep_install
    // downloads the tarball with Net::FTP unless --USE_HTTPS_PROTO is given, so on any
    // network that blocks passive FTP it dies with "Timeout at Net/FTP.pm" partway
    // through a ~27 GB transfer. --USE_HTTPS_PROTO is not a way out either: it requires
    // HTML::TableExtract, which this container does not ship, and vep_install responds
    // by printing "Skipping cache installation" and exiting 0 with no cache at all.
    //
    // Nothing is lost by skipping it. indexed_vep_cache tarballs are already the
    // tabix-converted form (all_vars.gz + .csi per chromosome), which is what the old
    // --CONVERT flag was there to produce, and they unpack straight into the
    // <species>/<version>_<assembly> layout VEP expects under --dir_cache.
    def filename = "${species}_vep_${cache_version}_${assembly}.tar.gz"
    def base_url = "https://ftp.ensembl.org/pub/release-${cache_version}/variation/indexed_vep_cache"
    """
    mkdir -p ${prefix}

    if [ "${preflight_check}" = "true" ]; then
        # CHECKSUMS is the release manifest. Fetching it up front turns a bad
        # species/assembly/version combination into a failure in seconds instead of one
        # discovered hours into the download, and it carries the checksum used below.
        wget -q -T 60 -O CHECKSUMS ${base_url}/CHECKSUMS

        awk -v f=${filename} '\$3 == f { found = 1 } END { exit !found }' CHECKSUMS || {
            echo "ERROR: ${filename} is not listed in ${base_url}/CHECKSUMS" >&2
            echo "Check --vep_species (${species}), --vep_genome (${assembly}) and --vep_cache_version (${cache_version})." >&2
            exit 1
        }
        echo "Pre-flight OK: ${filename} found in CHECKSUMS"
    fi

    # busybox wget has no --tries, so the retry loop is here. -c resumes the partial
    # file rather than restarting tens of GB from zero.
    downloaded=false
    for attempt in 1 2 3 4 5; do
        if wget -c -T 60 ${base_url}/${filename}; then
            downloaded=true
            break
        fi
        echo "Download attempt \$attempt for ${filename} failed, retrying" >&2
        sleep \$(( attempt * 30 ))
    done

    if [ "\$downloaded" != "true" ]; then
        echo "ERROR: could not download ${base_url}/${filename} after 5 attempts" >&2
        exit 1
    fi

    if [ "${preflight_check}" = "true" ]; then
        # Ensembl publishes BSD sum(1) checksums, "<checksum> <blocks> <filename>".
        # busybox sum omits the filename, so compare the first two fields, and compare
        # them as numbers because the checksum is zero-padded on one side only.
        expected=\$(awk -v f=${filename} '\$3 == f { print \$1 + 0, \$2 + 0 }' CHECKSUMS)
        observed=\$(sum ${filename} | awk '{ print \$1 + 0, \$2 + 0 }')

        if [ "\$expected" != "\$observed" ]; then
            echo "ERROR: checksum mismatch for ${filename} - the download is corrupt or truncated" >&2
            echo "  expected (CHECKSUMS): \$expected" >&2
            echo "  observed:             \$observed" >&2
            exit 1
        fi
        echo "Checksum OK: ${filename}"
    fi

    tar -xzf ${filename} -C ${prefix}
    rm -f ${filename} CHECKSUMS
    """

    stub:
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    mkdir -p ${prefix}/${species}/${cache_version}_${assembly}
    """
}
