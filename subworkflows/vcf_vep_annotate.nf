/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT ANNOTATION Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Annotates VCF files using Ensembl VEP.

    Cache resolution follows nf-core/sarek (subworkflows/nf-core/utils_annotation_cache):

      --download_cache      pull the cache with vep_install. One-off; publish it with
                            --outdir_cache and feed that path back as --vep_cache.
      --vep_cache <dir>     use an existing cache. This is the production path.

    download_cache wins when both are set, matching sarek's behaviour.

    VEP resolves a cache as <dir_cache>/<species>[_merged|_refseq]/<version>_<assembly>,
    so --vep_cache must point at the ROOT, not at the species or version directory. The
    layout is validated up front because vep's own failure mode is to silently fall back
    to a database connection and then die minutes later.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ENSEMBLVEP_DOWNLOAD } from '../modules/annotation/ensemblvep/download/main'
include { ENSEMBLVEP_VEP } from '../modules/annotation/ensemblvep/vep/main'
include { FASTVEP_ANNOTATE } from '../modules/annotation/fastvep/annotate/main'

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

    if (download_cache) {
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
    else if (vep_cache) {
        ch_vep_cache = channel.value(
            resolveVepCache(vep_cache, vep_genome, vep_species, vep_cache_version, params.vep_custom_args)
        )
    }
    else {
        error(
            "VEP annotation is enabled (--vep_mode ${params.vep_mode}) but no cache was given.\n" +
                "Pass --vep_cache <dir> to use an existing cache, or --download_cache to fetch one."
        )
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

    // FASTVEP is not implemented yet
    // FASTVEP_ANNOTATE()

    emit:
    vcf = ENSEMBLVEP_VEP.out.vcf // [meta, *_VEP.ann.vcf.gz]
    tab = ENSEMBLVEP_VEP.out.tab // [meta, *_VEP.ann.tab.gz]
    json = ENSEMBLVEP_VEP.out.json // [meta, *_VEP.ann.json.gz]
    reports = ENSEMBLVEP_VEP.out.report // *.summary.html
}

// Validate a user-supplied cache root and return the path VEP should be given as
// --dir_cache. Mirrors nf-core/sarek's UTILS_ANNOTATION_CACHE.
def resolveVepCache(cache, genome, species, cache_version, custom_args) {
    // Cloud-hosted caches (the annotation-cache buckets) nest the species directory
    // one level deeper, under <version>_<genome>/. Local caches do not.
    def cache_str = cache.toString()
    def is_cloud = ['s3://', 'gs://', 'az://'].any { cache_str.startsWith(it) }
    def cache_key = is_cloud ? "${cache_version}_${genome}/" : ''

    // --merged and --refseq caches live in differently named species directories.
    def args = custom_args ?: ''
    def species_suffix = args.contains('--merged') ? '_merged' : args.contains('--refseq') ? '_refseq' : ''

    def expected = "${cache_key}${species}${species_suffix}/${cache_version}_${genome}"
    def full = file("${cache_str}/${expected}", type: 'dir')

    if (!full.exists()) {
        error(
            "VEP cache not found: ${full}\n" +
                "--vep_cache must point at the cache ROOT, which VEP expects to contain\n" +
                "    ${expected}\n" +
                "Check --vep_cache_version (${cache_version}), --vep_genome (${genome}) and " +
                "--vep_species (${species}), or run with --download_cache to fetch it."
        )
    }

    return file("${cache_str}/${cache_key}", type: 'dir', checkIfExists: true)
}
