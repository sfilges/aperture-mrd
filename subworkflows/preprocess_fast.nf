/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAST PREPROCESS_READS Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Encapsulates the full read preprocessing chain:
      FastQC → fastp → alignment (BWA-MEM3/minibwa) → MarkDuplicates → indexing → QC

    Uses faster aligners and skips base quality score recalibration in comparison to 
    the 'standard' GATK workflow.

    Sample-type agnostic: tumor, normal, and plasma samples all pass through
    the same subworkflow. Branching by meta.status happens in main.nf.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW } from '../modules/fastqc/main'
include { FASTP } from '../modules/fastp/main'
include { PREPROCESS_ALIGN } from '../subworkflows/preprocess_align'
include { SAMTOOLS_MERGE } from '../modules/samtools/merge/main'
include { SAMTOOLS_MARKDUP } from '../modules/samtools/markdup/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from '../modules/samtools/index/main'
include { SAMTOOLS_STATS } from '../modules/samtools/stats/main'
include { RIKER_MULTI } from '../modules/fulcrum/riker/main'


workflow PREPROCESS_FAST {
    take:
    ch_reads // [meta, [fastq1, fastq2]]
    ch_fasta // [meta, fasta]
    ch_fasta_fai // [meta, fai]
    ch_index

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
            [groupKey(meta, meta.n_fastq), cram]
        }
        .groupTuple()

    //
    // Mark Duplicates
    //
    // Fastmode uses samtools instead of gatk, but SAMTOOLS_MARKDUP requires alignment files to be merged;
    // GATK4_MARKDUPLICATES can take a list of files
    SAMTOOLS_MERGE(
        cram_grouped,
        ch_fasta,
        ch_fasta_fai,
    )

    SAMTOOLS_MARKDUP(
        SAMTOOLS_MERGE.out.cram,
        ch_fasta,
        ch_fasta_fai,
    )

    //
    // Base Quality Score Recalibration
    //
    // Skip in fast mode

    //
    // Index recalibrated CRAM
    //
    INDEX_CRAM(
        SAMTOOLS_MARKDUP.out.cram
    )

    // Join CRAM with index → [meta, cram, crai]
    ch_cram = SAMTOOLS_MARKDUP.out.cram
        .join(
            INDEX_CRAM.out.crai,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, cram, crai -> [meta, cram, crai] }

    //
    // Alignment QC — riker
    //
    RIKER_MULTI(
        ch_cram,
        ch_fasta,
        ch_fasta_fai,
    )
    // TODO: Add reports to output channel

    emit:
    cram = ch_cram // [meta, cram, crai] — duplicate_marked, indexed
    reports = ch_reports // channel of QC files for MultiQC
}
