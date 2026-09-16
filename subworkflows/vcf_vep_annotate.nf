/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT ANNOTATION Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Annotates VCF files using Ensembl VEP.

    The cache is resolved in this order:

      --download_cache      force a fresh download, even when a usable cache is already
                            present.
      --vep_cache <dir>     use this cache. Must resolve, or the run fails: an explicit
                            path that does not validate is a typo, not a reason to start
                            a 27 GB download.
      (neither)             reuse the cache an earlier run downloaded into
                            <outdir_cache|outdir/cache>/vep_cache, and download it if
                            there is nothing to reuse.

    VEP resolves a cache as <dir_cache>/<species>[_merged|_refseq]/<version>_<assembly>,
    so --vep_cache must point at the ROOT, not at the species or version directory. The
    layout is validated up front because vep's own failure mode is to silently fall back
    to a database connection and then die minutes later.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ENSEMBLVEP_DOWNLOAD } from '../modules/ensemblvep/download/main'
include { ENSEMBLVEP_VEP      } from '../modules/ensemblvep/vep/main'

workflow VCF_VEP_ANNOTATE {
    take:
    vcf_for_vep // [meta, vcf]
    vep_genome // params.vep_genome
    vep_species // params.vep_species
    vep_cache_version // params.vep_cache_version
    vep_cache // params.vep_cache
    download_cache // params.download_cache
    fasta // [meta2, fasta]

    main:

    // ──────────────────────────────────────────────────────────────────────
    // Resolve the cache
    // ──────────────────────────────────────────────────────────────────────

    // Where an auto-downloaded cache ends up. Tracks the ENSEMBLVEP_DOWNLOAD publishDir
    // in conf/modules.config plus the module's own 'vep_cache' prefix.
    def published_cache = "${params.outdir_cache ?: "${params.outdir}/cache"}/vep_cache"

    def resolved_cache = null
    if (!download_cache) {
        if (vep_cache) {
            resolved_cache = resolveVepCache(vep_cache, vep_genome, vep_species, vep_cache_version, params.vep_custom_args)
        }
        else if (vepCacheExists(published_cache, vep_genome, vep_species, vep_cache_version, params.vep_custom_args)) {
            log.info("VEP: reusing the cache already downloaded to ${published_cache}")
            resolved_cache = resolveVepCache(published_cache, vep_genome, vep_species, vep_cache_version, params.vep_custom_args)
        }
    }

    if (resolved_cache) {
        ch_vep_cache = channel.value(resolved_cache)
    }
    else {
        ch_vep_info = channel.of(
            [
                ["id": "${vep_cache_version}_${vep_genome}"],
                vep_genome,
                vep_species,
                vep_cache_version,
            ]
        )

        ENSEMBLVEP_DOWNLOAD(ch_vep_info, params.vep_cache_preflight_check)

        // ENSEMBLVEP_VEP takes a bare `path cache`, so drop the meta and make it a
        // value channel — otherwise the single cache is consumed by the first VCF
        // and every later sample stalls.
        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map { _meta, cache -> cache }.first()
    }

    // ──────────────────────────────────────────────────────────────────────
    // Annotate
    // ──────────────────────────────────────────────────────────────────────

    // Method is controlled by modules.config via ext.when, only one method is used
    ENSEMBLVEP_VEP(
        vcf_for_vep,
        vep_genome,
        vep_species,
        vep_cache_version,
        ch_vep_cache,
        fasta,
    )

    emit:
    vcf = ENSEMBLVEP_VEP.out.vcf // [meta, *_VEP.ann.vcf.gz]
    tab = ENSEMBLVEP_VEP.out.tab // [meta, *_VEP.ann.tab.gz]
    json = ENSEMBLVEP_VEP.out.json // [meta, *_VEP.ann.json.gz]
    reports = ENSEMBLVEP_VEP.out.report // *.summary.html
}

// Where the cache data sits under a cache root, and what VEP should be handed as
// --dir_cache. Mirrors nf-core/sarek's UTILS_ANNOTATION_CACHE.
def vepCacheLayout(cache, genome, species, cache_version, custom_args) {
    // Cloud-hosted caches (the annotation-cache buckets) nest the species directory
    // one level deeper, under <version>_<genome>/. Local caches do not.
    def cache_str = cache.toString()
    def is_cloud = ['s3://', 'gs://', 'az://'].any { cache_str.startsWith(it) }
    def cache_key = is_cloud ? "${cache_version}_${genome}/" : ''

    // --merged and --refseq caches live in differently named species directories.
    def args = custom_args ?: ''
    def species_suffix = args.contains('--merged') ? '_merged' : args.contains('--refseq') ? '_refseq' : ''

    return [
        dir_cache: "${cache_str}/${cache_key}",
        expected: "${cache_key}${species}${species_suffix}/${cache_version}_${genome}",
    ]
}

// True when `cache` is a usable cache root. Used to probe for a cache an earlier run
// downloaded, where a miss means "download it" rather than "fail".
def vepCacheExists(cache, genome, species, cache_version, custom_args) {
    def layout = vepCacheLayout(cache, genome, species, cache_version, custom_args)
    return file("${cache}/${layout.expected}", type: 'dir').exists()
}

// Validate a cache root and return the path VEP should be given as --dir_cache.
// Unlike the probe above this fails hard, so a mistyped --vep_cache surfaces as an
// error instead of silently kicking off a download.
def resolveVepCache(cache, genome, species, cache_version, custom_args) {
    def layout = vepCacheLayout(cache, genome, species, cache_version, custom_args)
    def full = file("${cache}/${layout.expected}", type: 'dir')

    if (!full.exists()) {
        error(
            "VEP cache not found: ${full}\n" +
                "--vep_cache must point at the cache ROOT, which VEP expects to contain\n" +
                "    ${layout.expected}\n" +
                "Check --vep_cache_version (${cache_version}), --vep_genome (${genome}) and " +
                "--vep_species (${species}), or run with --download_cache to fetch it."
        )
    }

    return file(layout.dir_cache, type: 'dir', checkIfExists: true)
}
