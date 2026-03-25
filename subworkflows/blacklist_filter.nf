/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLACKLIST_FILTER Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-stage filtering of the somatic variant compendium:
      1. BED blacklist exclusion — removes variants in ENCODE blacklist,
         centromeric, and simple repeat regions
      2. gnomAD common variant exclusion — removes variants present in
         gnomAD (af-only-gnomad.hg38.vcf.gz) via bcftools isec

    VEP annotation is intentionally deferred to a separate annotation workflow
    to reduce computational cost during compendium building.

    Output: final patient SNV compendium VCF for MRDetect.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_CONCAT } from '../modules/bedtools/concat/main'
include { BEDTOOLS_SORT } from '../modules/bedtools/sort/main'
include { BEDTOOLS_MERGE } from '../modules/bedtools/merge/main'
include { BCFTOOLS_VIEW } from '../modules/bcftools/view/main'
include { BCFTOOLS_ISEC as GNOMAD_ISEC } from '../modules/bcftools/isec/main'

workflow BLACKLIST_FILTER {
    take:
    compendium_vcf // [meta, vcf.gz]
    compendium_tbi // [meta, tbi]
    gnomad_vcf // path: af-only-gnomad.hg38.vcf.gz (collected)
    gnomad_tbi // path: af-only-gnomad.hg38.vcf.gz.tbi (collected)
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
        BEDTOOLS_MERGE.out.bed.map { _meta, bed -> bed },
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    // ──────────────────────────────────────────────────────────────────────
    // Stage 2 — gnomAD common variant exclusion
    // Uses bcftools isec --complement to remove any variant present in
    // the gnomAD af-only VCF (all population variants, regardless of AF)
    // ──────────────────────────────────────────────────────────────────────

    // Build isec input: [meta, [compendium.vcf.gz, gnomad.vcf.gz], [compendium.tbi, gnomad.tbi]]
    ch_gnomad_isec_input = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .combine(gnomad_vcf)
        .combine(gnomad_tbi)
        .map { meta, vcf, tbi, gvcf, gtbi ->
            [meta, [vcf, gvcf], [tbi, gtbi]]
        }

    GNOMAD_ISEC(ch_gnomad_isec_input, "--complement")
    ch_versions = ch_versions.mix(GNOMAD_ISEC.out.versions)

    emit:
    vcf = GNOMAD_ISEC.out.vcf // [meta, compendium.filtered.vcf.gz]
    tbi = GNOMAD_ISEC.out.tbi // [meta, compendium.filtered.vcf.gz.tbi]
    versions = ch_versions
}
