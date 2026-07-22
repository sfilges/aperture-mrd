/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPROCESS_READS Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Encapsulates the full read preprocessing chain:
      FastQC → fastp → alignment → MarkDuplicates → BQSR → indexing → QC

    Sample-type agnostic: tumor, normal, and plasma samples all pass through
    the same subworkflow. Branching by meta.status happens in main.nf.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW } from '../modules/fastqc/main'
include { FASTP } from '../modules/fastp/main'
include { PREPROCESS_ALIGN } from '../subworkflows/preprocess_align'
include { GATK4_MARKDUPLICATES } from '../modules/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from '../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from '../modules/samtools/index/main'
include { SAMTOOLS_STATS } from '../modules/samtools/stats/main'
include { MOSDEPTH } from '../modules/mosdepth/main'
include { RIKER_MULTI } from '../modules/fulcrum/riker/main'

workflow PREPROCESS_GATK {
    take:
    ch_reads // [meta, [fastq1, fastq2]]
    ch_fasta // [meta, fasta]
    ch_fasta_fai // [meta, fai]
    ch_index // [meta, dir]
    dict // path(dict)
    known_sites_indels // [path(vcf), ...]
    known_sites_indels_tbi // [path(tbi), ...]

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
        params.split_fastq,
        params.trim_nextseq,
        params.length_required,
        params.save_fastqs,
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
    // Alignment — select aligner based on params.aligner
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
            [groupKey(meta, meta.n_fastq), cram]
        }
        .groupTuple()

    //
    // Mark Duplicates
    //
    GATK4_MARKDUPLICATES(
        cram_grouped,
        ch_fasta,
        ch_fasta_fai,
    )

    //
    // Base Quality Score Recalibration
    //
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

    //
    // Index recalibrated CRAM
    //
    INDEX_CRAM(
        GATK4_APPLYBQSR.out.cram
    )

    // Join CRAM with index → [meta, cram, crai]
    ch_cram = GATK4_APPLYBQSR.out.cram
        .join(
            INDEX_CRAM.out.crai,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, cram, crai -> [meta, cram, crai] }

    //
    // Alignment QC — riker
    //
    RIKER_MULTI(ch_cram, ch_fasta, ch_fasta_fai)
    // TODO: Add reports to output channel

    emit:
    cram = ch_cram // [meta, cram, crai] — recalibrated, indexed
    reports = ch_reports // channel of QC files for MultiQC
}
