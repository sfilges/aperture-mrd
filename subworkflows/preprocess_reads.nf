/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPROCESS_READS Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Encapsulates the full read preprocessing chain:
      FastQC → fastp → alignment → MarkDuplicates → [BQSR] → indexing → QC

    params.preprocessing selects the dedup/recalibration strategy:
      'gatk' : GATK4 MarkDuplicates → BaseRecalibrator → ApplyBQSR
      'fast' : samtools merge + markdup, skips BQSR — tuned for bwamem3/minibwa

    Sample-type agnostic: tumor, normal, and plasma samples all pass through
    the same subworkflow. Branching by meta.status happens in main.nf.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW } from '../modules/fastqc/main'
include { FASTP } from '../modules/fastp/main'
include { PREPROCESS_ALIGN } from '../subworkflows/preprocess_align'
include { SAMTOOLS_MERGE } from '../modules/samtools/merge/main'
include { SAMTOOLS_MARKDUP } from '../modules/samtools/markdup/main'
include { GATK4_MARKDUPLICATES } from '../modules/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from '../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from '../modules/samtools/index/main'
include { SAMTOOLS_STATS } from '../modules/samtools/stats/main'
include { RIKER_MULTI } from '../modules/fulcrum/riker/multi/main'

workflow PREPROCESS_READS {
    take:
    ch_reads // [meta, [fastq1, fastq2]]
    ch_fasta // [meta, fasta]
    ch_fasta_fai // [meta, fai]
    ch_index // [meta, dir]
    mode // 'wgs' or 'wes'
    baits_in // channel: bait intervals — wes only
    targets_in // channel: target intervals — wes only
    dict // path(dict) — gatk only
    known_sites_indels // [path(vcf), ...] — gatk only
    known_sites_indels_tbi // [path(tbi), ...] — gatk only

    main:
    ch_reports = channel.empty()

    //
    // FastQC on raw reads
    //
    FASTQC_RAW(ch_reads)
    ch_reports = ch_reports.mix(FASTQC_RAW.out.zip.collect { it -> it[1] }.ifEmpty([]))

    //
    // Trimming with fastp
    //
    FASTP(
        ch_reads,
        params.merge_fastqs,
    )
    ch_reports = ch_reports.mix(FASTP.out.json.collect { it -> it[1] }.ifEmpty([]))

    // Handle split FASTQs: collate paired reads and transpose for per-chunk alignment
    if (params.split_fastq) {
        ch_trim_reads = FASTP.out.reads
            .map { meta, reads ->
                def read_files = reads.sort(false) { a, b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                [meta + [n_fastq: read_files.size()], read_files]
            }
            .transpose()
    }
    else {
        ch_trim_reads = FASTP.out.reads
    }

    if (params.merge_fastqs) {
        ch_trim_reads = ch_trim_reads.mix(FASTP.out.reads_merged)
    }

    //
    // Alignment — select aligner based on params.aligner + ext.when
    //
    PREPROCESS_ALIGN(
        ch_trim_reads,
        ch_fasta,
        ch_fasta_fai,
        ch_index,
    )

    //
    // Group split CRAMs back by sample before deduplication
    //
    cram_grouped = PREPROCESS_ALIGN.out.crams
        .map { meta, cram ->
            // n_fastq is only set when split_fastq is enabled; default to 1 chunk otherwise
            [groupKey(meta, meta.n_fastq ?: 1), cram]
        }
        .groupTuple()

    //
    // Mark Duplicates [+ BQSR for gatk]
    //
    if (params.preprocessing == 'gatk') {
        GATK4_MARKDUPLICATES(
            cram_grouped,
            ch_fasta,
            ch_fasta_fai,
        )

        GATK4_BASERECALIBRATOR(
            GATK4_MARKDUPLICATES.out.cram,
            ch_fasta,
            ch_fasta_fai,
            dict,
            known_sites_indels,
            known_sites_indels_tbi,
        )
        ch_reports = ch_reports.mix(GATK4_BASERECALIBRATOR.out.table.collect { _meta, table -> [table] })

        // Join CRAM with recalibration table for ApplyBQSR
        cram_applybqsr = GATK4_MARKDUPLICATES.out.cram.join(GATK4_BASERECALIBRATOR.out.table, failOnDuplicate: true, failOnMismatch: true)

        GATK4_APPLYBQSR(
            cram_applybqsr,
            ch_fasta,
            ch_fasta_fai,
            dict,
        )

        ch_cram_dedup = GATK4_APPLYBQSR.out.cram
    }
    else {
        // Fast mode uses samtools instead of gatk, and SAMTOOLS_MARKDUP takes a single
        // input file (GATK4_MARKDUPLICATES can take a list directly). Merge is only needed
        // when a sample was split into multiple chunks; single-chunk samples (e.g. when
        // split_fastq is off) skip the no-op merge and go straight to markdup. Skips BQSR.
        cram_grouped
            .branch { _meta, crams ->
                single: crams.size() == 1
                multi: true
            }
            .set { ch_cram_to_merge }

        SAMTOOLS_MERGE(
            ch_cram_to_merge.multi,
            ch_fasta,
            ch_fasta_fai,
        )

        ch_cram_to_markdup = SAMTOOLS_MERGE.out.cram
            .mix(ch_cram_to_merge.single.map { meta, crams -> [meta, crams[0]] })

        SAMTOOLS_MARKDUP(
            ch_cram_to_markdup,
            ch_fasta,
            ch_fasta_fai,
        )

        ch_cram_dedup = SAMTOOLS_MARKDUP.out.cram
    }

    //
    // Index CRAM
    //
    INDEX_CRAM(
        ch_cram_dedup
    )

    // Join CRAM with index → [meta, cram, crai]
    ch_cram = ch_cram_dedup
        .join(
            INDEX_CRAM.out.crai,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, cram, crai -> [meta, cram, crai] }

    //
    // Alignment QC — riker
    //
    ch_hybcap = (mode == 'wes')
        ? baits_in.combine(targets_in).map { baits, targets -> [["id": "wes_baits_targets"], baits, targets] }
        : channel.value([[:], [], []])

    RIKER_MULTI(
        ch_cram,
        ch_fasta,
        ch_fasta_fai,
        ch_hybcap,
    )
    ch_reports = ch_reports.mix(
        RIKER_MULTI.out.alignment_metrics.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.base_dist.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.mean_qual.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.qual_dist.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.isize_metrics.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.isize_histogram.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.gcbias_detail.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.gcbias_summary.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.error_indel.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.error_mismatch.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.error_overlap.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.wgs_coverage.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.wgs_metrics.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.hybcap_metrics.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.hybcap_per_target.collect { it -> it[1] }.ifEmpty([]),
        RIKER_MULTI.out.hybcap_per_base.collect { it -> it[1] }.ifEmpty([]),
    )

    emit:
    cram = ch_cram // [meta, cram, crai] — duplicate_marked[, recalibrated], indexed
    reports = ch_reports // channel of QC files for MultiQC
}
