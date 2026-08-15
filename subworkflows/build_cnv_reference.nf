/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BUILD_CNV_REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Build the assay-level copy number reference bundle from a set of process-
    matched normal samples. Run once per assay, not per cohort:

        nextflow run . --step build_cnv_reference --input normals.csv

    Produces
      <name>.reference.cnn      CNVkit pooled reference  -> --cnvkit_reference
      normalDB_<name>_<genome>.rds   PureCN coverage database -> --purecn_normaldb
      mapping_bias_<name>_<genome>.rds  PureCN mapping bias   -> --purecn_mapping_bias

    A pooled reference is the recommended input for CNVkit: it captures the
    systematic coverage biases of the assay far better than a single matched
    normal, and much better than the flat fallback.

    CNVkit builds the panel from all normals in one batch call with no tumor
    positional argument, which also writes the per-normal coverage that PureCN's
    NormalDB.R consumes directly — no separate PureCN coverage step is needed,
    since NormalDB.R accepts CNVkit format natively.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNVKIT_ACCESS    } from '../modules/cnvkit/access/main'
include { CNVKIT_BATCH     } from '../modules/cnvkit/batch/main'
include { PURECN_NORMALDB  } from '../modules/purecn/normaldb/main'

workflow BUILD_CNV_REFERENCE {
    take:
    ch_normal_bams // [meta, bam, bai] every normal in the panel
    ch_fasta // [meta, fasta]
    ch_fai // [meta, fai]
    ch_targets // [meta, baits.bed] — probes, unpadded
    ch_exclude // [bed, ...] regions to drop from the accessible genome

    main:

    def reference_name = params.cnv_reference_name

    //
    // Accessible regions — genome-level, built once
    //
    if (params.cnvkit_access) {
        ch_access = channel.value([["id": "access"], file(params.cnvkit_access, checkIfExists: true)])
    }
    else {
        CNVKIT_ACCESS(ch_fasta, ch_exclude)
        ch_access = CNVKIT_ACCESS.out.bed.first()
    }

    ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { meta, fasta, _meta2, fai -> [meta, fasta, fai] }
        .first()

    //
    // CNVkit — pooled reference from every normal at once
    //
    // Collapsing all normals under one synthetic meta is what puts CNVKIT_BATCH
    // into panel mode: no tumor positional argument, and every normal passed to
    // a single --normal.
    ch_pon_input = ch_normal_bams
        .map { _meta, bam, bai -> [bam, bai] }
        .collect(flat: false)
        .map { pairs ->
            [["id": reference_name], [], [], pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] }]
        }

    CNVKIT_BATCH(
        ch_pon_input,
        ch_fasta_fai,
        ch_targets,
        ch_access,
        channel.value([[:], []]),
        params.mode == 'wgs' ? 'wgs' : 'hybrid',
    )

    //
    // PureCN — coverage database and, given a panel VCF, the mapping bias
    //
    // Only the on-target coverage is passed: PureCN's CNVkit reader marks every
    // interval it is given as on-target, so including the antitarget files would
    // silently mislabel off-target bins.
    ch_normal_panel = params.purecn_normal_panel
        ? channel.value([
            file(params.purecn_normal_panel, checkIfExists: true),
            file("${params.purecn_normal_panel}.tbi", checkIfExists: true),
        ])
        : channel.value([[], []])

    ch_normaldb_input = CNVKIT_BATCH.out.target_cnn
        .combine(ch_normal_panel)
        .map { meta, coverages, vcf, tbi -> [meta, coverages, vcf, tbi] }

    PURECN_NORMALDB(
        ch_normaldb_input,
        params.purecn_genome,
        reference_name,
    )

    emit:
    reference_cnn = CNVKIT_BATCH.out.reference_cnn // [meta, .reference.cnn]
    target_cnn = CNVKIT_BATCH.out.target_cnn // [meta, [*.targetcoverage.cnn]]
    access_bed = ch_access // [meta, access.bed]
    normaldb = PURECN_NORMALDB.out.rds // [meta, normalDB_*.rds]
    mapping_bias = PURECN_NORMALDB.out.bias_rds // [meta, mapping_bias_*.rds]
    mapping_bias_bed = PURECN_NORMALDB.out.bias_bed // [meta, mapping_bias_hq_sites_*.bed]
    low_coverage_bed = PURECN_NORMALDB.out.low_cov_bed // [meta, low_coverage_targets_*.bed]
    interval_weights = PURECN_NORMALDB.out.png // [meta, interval_weights_*.png]
}
