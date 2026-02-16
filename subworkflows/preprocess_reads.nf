/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPROCESS_READS Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Encapsulates the full read preprocessing chain:
      FastQC → fastp → alignment (BWA/BWA-MEM2) → MarkDuplicates → BQSR → indexing → QC

    Sample-type agnostic: tumor, normal, and plasma samples all pass through
    the same subworkflow. Branching by meta.status happens in main.nf.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW     } from '../modules/fastqc/main'
include { FASTP                     } from '../modules/fastp/main'
include { BWA_INDEX                 } from '../modules/bwamem/index/main'
include { BWA_MEM                   } from '../modules/bwamem/mem/main'
include { BWAMEM2_INDEX             } from '../modules/bwamem2/index/main'
include { BWAMEM2_MEM              } from '../modules/bwamem2/mem/main'
include { GATK4_MARKDUPLICATES      } from '../modules/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR    } from '../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR           } from '../modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from '../modules/samtools/index/main'
include { SAMTOOLS_STATS            } from '../modules/samtools/stats/main'
include { MOSDEPTH                  } from '../modules/mosdepth/main'

workflow PREPROCESS_READS {

    take:
    ch_reads               // [meta, [fastq1, fastq2]]
    ch_fasta               // [meta, fasta]
    ch_fasta_fai           // [meta, fai]
    dict                   // path(dict)
    known_sites_indels     // [path(vcf), ...]
    known_sites_indels_tbi // [path(tbi), ...]

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    //
    // FastQC on raw reads
    //
    FASTQC_RAW(ch_reads)
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)
    ch_reports  = ch_reports.mix(FASTQC_RAW.out.zip.collect { it[1] }.ifEmpty([]))

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
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_reports  = ch_reports.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))

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
    if (params.aligner == 'bwa') {
        if (params.bwaindex) {
            ch_index = channel.value(
                [["id": "bwa_index"], file("${params.bwaindex}/*{,.amb,.ann,.bwt,.pac,.sa}")]
            )
        }
        else {
            BWA_INDEX(ch_fasta)
            ch_index = BWA_INDEX.out.index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }

        BWA_MEM(
            ch_trim_reads,
            ch_index,
        )

        ch_aligned_bams = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    }
    else if (params.aligner == 'bwamem2') {
        if (params.bwamem2index) {
            ch_index = channel.value(
                [["id": "bwamem2_index"], file("${params.bwamem2index}/*{,.amb,.ann,.bwt,.pac,.sa}")]
            )
        }
        else {
            BWAMEM2_INDEX(ch_fasta)
            ch_index = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }

        BWAMEM2_MEM(
            ch_trim_reads,
            ch_fasta,
            ch_index,
        )

        ch_aligned_bams = BWAMEM2_MEM.out.bam
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    }
    else {
        error("Unknown aligner: ${params.aligner}. Supported aligners are: bwa, bwamem2.")
    }

    //
    // Group split BAMs back by sample before deduplication
    //
    bam_grouped = ch_aligned_bams
        .map { meta, bam ->
            [groupKey(meta, meta.n_fastq), bam]
        }
        .groupTuple()

    //
    // Mark Duplicates
    //
    GATK4_MARKDUPLICATES(
        bam_grouped,
        ch_fasta,
        ch_fasta_fai,
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

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
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)
    ch_reports  = ch_reports.mix(GATK4_BASERECALIBRATOR.out.table.collect { meta, table -> [table] })

    // Join CRAM with recalibration table for ApplyBQSR
    cram_applybqsr = GATK4_MARKDUPLICATES.out.cram
        .join(GATK4_BASERECALIBRATOR.out.table, failOnDuplicate: true, failOnMismatch: true)

    GATK4_APPLYBQSR(
        cram_applybqsr,
        ch_fasta,
        ch_fasta_fai,
        dict,
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

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
    // Alignment QC — samtools stats + mosdepth
    //
    SAMTOOLS_STATS(ch_cram, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_reports  = ch_reports.mix(SAMTOOLS_STATS.out.stats.collect { meta, stats -> [stats] })

    MOSDEPTH(ch_cram, ch_fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_reports  = ch_reports.mix(MOSDEPTH.out.global_txt.collect { meta, global_txt -> [global_txt] })
    ch_reports  = ch_reports.mix(MOSDEPTH.out.regions_txt.collect { meta, regions_txt -> [regions_txt] })

    emit:
    cram     = ch_cram      // [meta, cram, crai] — recalibrated, indexed
    reports  = ch_reports    // channel of QC files for MultiQC
    versions = ch_versions   // channel of versions.yml
}
