/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC_CNV_CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Somatic copy number calling with CNVkit, optionally refined by PureCN.

    Handles tumor-normal and tumor-only in one flow. Which copy number reference
    is used is decided once, up front:

      params.cnvkit_reference set  -> that pooled reference, matched normals ignored
      matched normal present       -> paired reference built from that normal
      neither                      -> flat reference (equal coverage assumed)

    The flat reference is the weakest of the three and is a fallback, not a
    recommendation — build a panel with build_cnv_reference where possible.

    With PureCN enabled the CNVkit segments are re-called using its purity and
    ploidy estimate, which is materially better than fixed log2 thresholds.

    Inputs are BAM. CNVkit reads CRAM natively but is prohibitively slow doing
    so, so conversion happens upstream.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNVKIT_ACCESS     } from '../modules/cnvkit/access/main'
include { CNVKIT_BATCH      } from '../modules/cnvkit/batch/main'
include { CNVKIT_EXPORT_SEG } from '../modules/cnvkit/export/main'
include { CNVKIT_CALL       } from '../modules/cnvkit/call/main'
include { PURECN_RUN        } from '../modules/purecn/run/main'

workflow SOMATIC_CNV_CALLING {
    take:
    ch_bam_pairs // [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]; normal may be []
    ch_fasta // [meta, fasta]
    ch_fai // [meta, fai]
    ch_targets // [meta, baits.bed] — probes, unpadded
    ch_vcf // [meta, vcf, tbi] Mutect2 calls including germline sites
    ch_exclude // [bed, ...] regions to drop from the accessible genome

    main:

    // Setting --cnvkit_reference discards the matched normal for EVERY sample in
    // the run, including tumor/normal pairs, which are then handled exactly as
    // tumor-only. cnvkit cannot use both anyway (it exits 1 when -n, -t or -g is
    // passed alongside -r), so this is a per-run choice rather than a per-sample
    // one: to keep matched normals for the pairs, leave --cnvkit_reference unset.
    def use_reference = params.cnvkit_reference as boolean

    ch_reference = use_reference
        ? channel.value([["id": "cnvkit_reference"], file(params.cnvkit_reference, checkIfExists: true)])
        : channel.value([[:], []])

    //
    // Accessible regions — genome-level, built once. Not needed at all when
    // reusing a reference, since the bins are already fixed by the .cnn.
    //
    if (use_reference) {
        ch_access = channel.value([[:], []])
    }
    else if (params.cnvkit_access) {
        ch_access = channel.value([["id": "access"], file(params.cnvkit_access, checkIfExists: true)])
    }
    else {
        CNVKIT_ACCESS(ch_fasta, ch_exclude)
        ch_access = CNVKIT_ACCESS.out.bed.first()
    }

    //
    // CNVkit batch — one call covers all three reference modes
    //
    ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { meta, fasta, _meta2, fai -> [meta, fasta, fai] }
        .first()

    ch_batch_input = ch_bam_pairs.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
        use_reference
            ? [meta, tumor_bam, tumor_bai, [], []]
            : [meta, tumor_bam, tumor_bai, normal_bam ?: [], normal_bai ?: []]
    }

    CNVKIT_BATCH(
        ch_batch_input,
        ch_fasta_fai,
        ch_targets,
        ch_access,
        ch_reference,
        params.mode == 'wgs' ? 'wgs' : 'hybrid',
    )

    //
    // PureCN — purity and ploidy from the CNVkit output
    //
    // PureCN consumes .cnr as --tumor plus the DNAcopy-format .seg as
    // --seg-file, which is why it needs neither an interval file nor coverage
    // files of its own.
    //
    CNVKIT_EXPORT_SEG(CNVKIT_BATCH.out.cns)

    if (params.purecn) {
        ch_purecn_input = CNVKIT_BATCH.out.cnr
            .join(CNVKIT_EXPORT_SEG.out.seg, failOnDuplicate: true, failOnMismatch: true)
            .join(ch_vcf, failOnDuplicate: true, failOnMismatch: true)

        // The SNP blacklist matters only when neither a matched normal nor a
        // pool of normals is available. This input is assay-level rather than
        // per-sample, so a configured pooled reference is the only signal
        // available here; per-pair matched-normal state cannot be seen.
        ch_snp_blacklist = !use_reference && params.simple_repeats
            ? file(params.simple_repeats, checkIfExists: true)
            : []

        PURECN_RUN(
            ch_purecn_input,
            [[:], []],
            params.purecn_normaldb ? file(params.purecn_normaldb, checkIfExists: true) : [],
            params.purecn_mapping_bias ? file(params.purecn_mapping_bias, checkIfExists: true) : [],
            ch_snp_blacklist,
            params.purecn_genome,
        )

        // The curation CSV carries one row per sample with Purity and Ploidy
        ch_purity_ploidy = PURECN_RUN.out.csv
            .splitCsv(header: true, elem: 1)
            .map { meta, row -> [meta, row.Purity, row.Ploidy] }

        ch_call_input = CNVKIT_BATCH.out.cns.join(ch_purity_ploidy, failOnDuplicate: true, failOnMismatch: true)

        ch_purecn_csv = PURECN_RUN.out.csv
        ch_purecn_seg = PURECN_RUN.out.seg
    }
    else {
        // No purity estimate: CNVKIT_CALL falls back to fixed log2 thresholds
        ch_call_input = CNVKIT_BATCH.out.cns.map { meta, cns -> [meta, cns, null, null] }
        ch_purecn_csv = channel.empty()
        ch_purecn_seg = channel.empty()
    }

    CNVKIT_CALL(ch_call_input)

    emit:
    cnr = CNVKIT_BATCH.out.cnr // [meta, .cnr]  bin-level log2 ratios
    cns = CNVKIT_BATCH.out.cns // [meta, .cns]  segments
    called_cns = CNVKIT_CALL.out.cns // [meta, .called.cns] integer copy number
    bintest_cns = CNVKIT_BATCH.out.bintest_cns // [meta, .bintest.cns]
    reference_cnn = CNVKIT_BATCH.out.reference_cnn // [meta, .reference.cnn] when built
    seg = CNVKIT_EXPORT_SEG.out.seg // [meta, .seg]  DNAcopy format
    purecn_csv = ch_purecn_csv // [meta, .csv]  purity/ploidy
    purecn_seg = ch_purecn_seg // [meta, _dnacopy.seg]
}
