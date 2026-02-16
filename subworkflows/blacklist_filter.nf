/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLACKLIST_FILTER Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three-stage filtering of the somatic variant compendium:
      1. BED blacklist exclusion — removes variants in ENCODE blacklist,
         centromeric, and simple repeat regions
      2. VEP annotation — adds consequence, gene, gnomAD allele frequencies
      3. gnomAD AF filter — removes common variants (gnomADg_AF > 0.01)

    Output: final patient SNV compendium VCF for MRDetect.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_CONCAT } from '../modules/bedtools/concat/main'
include { BEDTOOLS_SORT } from '../modules/bedtools/sort/main'
include { BEDTOOLS_MERGE } from '../modules/bedtools/merge/main'
include { BCFTOOLS_VIEW } from '../modules/bcftools/view/main'
include { ENSEMBLVEP_VEP } from '../modules/ensemblvep/vep/main'
include { BCFTOOLS_FILTER as GNOMAD_FILTER } from '../modules/bcftools/filter/main'

workflow BLACKLIST_FILTER {
    take:
    compendium_vcf // [meta, vcf.gz]
    compendium_tbi // [meta, tbi]
    ch_fasta // [meta, fasta]
    ch_fasta_fai // [meta, fai]
    vep_genome // val: 'GRCh38'
    vep_species // val: 'homo_sapiens'
    vep_cache_version // val: '111'
    vep_cache // path: cache dir (or [])
    blacklists // [bed1, bed2, ...] collected

    main:
    ch_versions = channel.empty()

    // ──────────────────────────────────────────────────────────────────────
    // Stage 1 — BED blacklist exclusion
    // ──────────────────────────────────────────────────────────────────────

    // Concatenate all blacklist BEDs into one, then sort and merge
    ch_beds = channel.value([["id": "blacklist"], blacklists])

    BEDTOOLS_CONCAT(ch_beds)
    ch_versions = ch_versions.mix(BEDTOOLS_CONCAT.out.versions)

    BEDTOOLS_SORT(BEDTOOLS_CONCAT.out.bed, [])
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    BEDTOOLS_MERGE(BEDTOOLS_SORT.out.sorted)
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

    // Join VCF + TBI, then exclude blacklisted regions
    ch_view_input = compendium_vcf.join(compendium_tbi, failOnDuplicate: true, failOnMismatch: true)

    BCFTOOLS_VIEW(
        ch_view_input,
        BEDTOOLS_MERGE.out.bed.map { meta, bed -> bed },
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    // ──────────────────────────────────────────────────────────────────────
    // Stage 2 — VEP annotation (skipped if no cache provided)
    // ──────────────────────────────────────────────────────────────────────

    // Prepare VEP input: [meta, vcf, extra_files]
    ch_vep_input = BCFTOOLS_VIEW.out.vcf.map { meta, vcf -> [meta, vcf, []] }

    ENSEMBLVEP_VEP(
        ch_vep_input,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        ch_fasta,
        [],
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    // ──────────────────────────────────────────────────────────────────────
    // Stage 3 — gnomAD AF filter (expression via modules.config ext.args)
    // ──────────────────────────────────────────────────────────────────────

    GNOMAD_FILTER(ENSEMBLVEP_VEP.out.vcf)
    ch_versions = ch_versions.mix(GNOMAD_FILTER.out.versions)

    emit:
    vcf = GNOMAD_FILTER.out.vcf // [meta, compendium.filtered.vcf.gz]
    tbi = GNOMAD_FILTER.out.tbi // [meta, compendium.filtered.vcf.gz.tbi]
    versions = ch_versions
}
