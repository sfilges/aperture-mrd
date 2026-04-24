// Mutect2 and post-processing modules
include { GATK4_MUTECT2 as MUTECT2_SOMATIC } from '../modules/gatk4/mutect2/main'
include { GATK4_MERGEVCFS as MERGE_MUTECT2 } from '../modules/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS as MERGEMUTECTSTATS } from '../modules/gatk4/mergemutectstats/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../modules/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL } from '../modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR } from '../modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION } from '../modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS as FILTERMUTECTCALLS } from '../modules/gatk4/filtermutectcalls/main'

// Other variant callers
include { MANTA_SOMATIC } from '../modules/variant_callers/manta/main'
include { STRELKA_SOMATIC } from '../modules/variant_callers/strelka2/main'
include { LOFREQ_SOMATIC } from '../modules/variant_callers/lofreq/main'
include { MUSE_CALL } from '../modules/variant_callers/muse/call/main'

// Utility modules
include { TABIX_BGZIPTABIX } from '../modules/tabix/tabix_bgzip/main'
include { BEDTOOLS_SPLIT } from '../modules/bedtools/split/main'

workflow TN_SOMATIC_VARIANT_CALLING {
    take:
    cram_variant_calling_pair
    ch_fasta
    ch_fasta_fai
    dict
    germline_resource
    germline_resource_tbi
    versions
    dbsnp
    dbsnp_tbi
    _intervals_bed_all
    intervals_bed_gbz_tbi_all
    intervals_bed_split
    intervals_bed_bgz_tbi_split // [meta, intervals.bed, intervals.bed.tbi]

    main:

    // Wrap dict as tuple for modules that require [meta, dict] format
    ch_dict_meta = dict.map { it -> [["id": "dict"], it] }

    // =========================================================================
    // MUTECT2 — Scatter by intervals
    // =========================================================================

    // Scatter: combine each tumor-normal pair with each interval BED for parallel execution
    cram_variant_calling_pair_interval_bed_split = cram_variant_calling_pair.combine(intervals_bed_split.map { _meta, bed -> [bed] })

    // Channel: [[meta], normal_cram, normal_crai, tumor_cram, tumor_crai, intervals.bed]

    MUTECT2_SOMATIC(
        cram_variant_calling_pair_interval_bed_split,
        ch_fasta,
        ch_fasta_fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        [],
        [],
    )

    // =========================================================================
    // MUTECT2 POST-PROCESSING — Gather scattered outputs and filter
    // =========================================================================

    // Gather: group scattered outputs back by sample pair
    MERGE_MUTECT2(
        MUTECT2_SOMATIC.out.vcf.groupTuple(),
        ch_dict_meta,
    )

    MERGEMUTECTSTATS(
        MUTECT2_SOMATIC.out.stats.groupTuple()
    )

    // Learn read orientation model from f1r2 counts (used to filter orientation bias artifacts)
    LEARNREADORIENTATIONMODEL(
        MUTECT2_SOMATIC.out.f1r2.groupTuple()
    )

    // GetPileupSummaries for tumor and normal — needed for contamination estimation
    // Extract tumor and normal CRAMs separately from the pair channel
    ch_tumor_pileup_input = cram_variant_calling_pair.map { meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai ->
        [meta, tumor_cram, tumor_crai, []]
    }

    ch_normal_pileup_input = cram_variant_calling_pair.map { meta, normal_cram, normal_crai, _tumor_cram, _tumor_crai ->
        [meta, normal_cram, normal_crai, []]
    }

    GETPILEUPSUMMARIES_TUMOR(
        ch_tumor_pileup_input,
        ch_fasta,
        ch_fasta_fai,
        ch_dict_meta,
        germline_resource,
        germline_resource_tbi,
    )

    GETPILEUPSUMMARIES_NORMAL(
        ch_normal_pileup_input,
        ch_fasta,
        ch_fasta_fai,
        ch_dict_meta,
        germline_resource,
        germline_resource_tbi,
    )

    // Calculate contamination: tumor pileup with matched normal pileup
    ch_contamination_input = GETPILEUPSUMMARIES_TUMOR.out.table.join(GETPILEUPSUMMARIES_NORMAL.out.table, failOnDuplicate: true, failOnMismatch: true)

    CALCULATECONTAMINATION(ch_contamination_input)

    // Assemble all inputs for FilterMutectCalls
    // Required input: [meta, vcf, tbi, stats, orientationbias, segmentation, contamination_table, contamination_estimate]
    ch_filter_input = MERGE_MUTECT2.out.vcf
        .join(MERGE_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(MERGEMUTECTSTATS.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)
        .join(CALCULATECONTAMINATION.out.segmentation, failOnDuplicate: true, failOnMismatch: true)
        .join(CALCULATECONTAMINATION.out.contamination, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi, stats, artifactprior, segmentation, contamination ->
            [meta, vcf, tbi, stats, artifactprior, segmentation, contamination, 0]
        }

    FILTERMUTECTCALLS(
        ch_filter_input,
        ch_fasta,
        ch_fasta_fai,
        ch_dict_meta,
    )

    // =========================================================================
    // MANTA and STRELKA
    // =========================================================================

    // Run Manta for somatic SV variant calling and candidate small indels for Strelka
    // Provide all intervals in one file, bg-zipped and tabix indexed
    MANTA_SOMATIC(
        cram_variant_calling_pair,
        ch_fasta,
        ch_fasta_fai,
        intervals_bed_gbz_tbi_all,
        [],
    )

    // Collect the output candidate small indels VCF and TBI files from Manta for Strelka
    cram_strelka = cram_variant_calling_pair
        .join(MANTA_SOMATIC.out.candidate_small_indels_vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi, failOnDuplicate: true, failOnMismatch: true)
        .combine(intervals_bed_bgz_tbi_split.map { _meta, bed, tbi -> [bed, tbi] })

    // Strelka calls the entire genome by default, however variant calling may be restricted to an arbitrary
    // subset of the genome by providing a region file in BED format with the --callRegions configuration option.
    // The BED file must be bgzip-compressed and tabix-indexed, and only one such BED file may be specified
    STRELKA_SOMATIC(
        cram_strelka,
        ch_fasta,
        ch_fasta_fai,
    )

    // =========================================================================
    // LOFREQ
    // =========================================================================

    // TODO: Convert CRAM to BAM for LOFREQ_SOMATIC and MUSE_SOMATIC
    LOFREQ_SOMATIC(
        cram_variant_calling_pair,
        ch_fasta,
        ch_fasta_fai,
        dbsnp,
        dbsnp_tbi,
    )

    //
    // MUSE
    //

    //MUSE_SOMATIC(
    //    bam_variant_calling_pair,
    //    ch_fasta
    //)

    // =========================================================================
    // COLLECT OUTPUTS
    // =========================================================================

    versions = versions.mix(MUTECT2_SOMATIC.out.versions)
    versions = versions.mix(MERGE_MUTECT2.out.versions)
    versions = versions.mix(MERGEMUTECTSTATS.out.versions)
    versions = versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    versions = versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions)
    versions = versions.mix(CALCULATECONTAMINATION.out.versions)
    versions = versions.mix(FILTERMUTECTCALLS.out.versions)
    versions = versions.mix(MANTA_SOMATIC.out.versions)
    versions = versions.mix(STRELKA_SOMATIC.out.versions)
    versions = versions.mix(LOFREQ_SOMATIC.out.versions)

    emit:
    mutect2_vcf = FILTERMUTECTCALLS.out.vcf // [meta, filtered.vcf.gz]
    mutect2_tbi = FILTERMUTECTCALLS.out.tbi // [meta, filtered.vcf.gz.tbi]
    mutect2_stats = FILTERMUTECTCALLS.out.stats // [meta, filteringStats.tsv]
    strelka_snvs_vcf = STRELKA_SOMATIC.out.vcf_snvs // [meta, somatic_snvs.vcf.gz]
    strelka_snvs_tbi = STRELKA_SOMATIC.out.vcf_snvs_tbi // [meta, somatic_snvs.vcf.gz.tbi]
    strelka_indels_vcf = STRELKA_SOMATIC.out.vcf_indels // [meta, somatic_indels.vcf.gz]
    strelka_indels_tbi = STRELKA_SOMATIC.out.vcf_indels_tbi // [meta, somatic_indels.vcf.gz.tbi]
    manta_sv_vcf = MANTA_SOMATIC.out.somatic_sv_vcf // [meta, somaticSV.vcf.gz]
    lofreq_vcf = LOFREQ_SOMATIC.out.vcf // [meta, *.vcf.gz]
    versions = versions
}
