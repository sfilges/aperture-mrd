/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMI Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fulcrum fgumi implementation (Rust-based, supersedes previous fgbio).

    Bypasses standard preprocxessign workflows ('fast', 'gatk') to produce 
    a bam file of aligned and sorted UMI-corrected reads from fastq inputs.

        - Uses bwa-mem3 as the aligner
        - Output cram

    // TODO: Or should the subprocess only output UMI-corrected reads as fastq
    // for further processing according to the standard workflows above?

    Source: https://github.com/fulcrumgenomics/fgumi
    License: MIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FULCRUM_FGUMI_EXTRACT } from '../modules/fulcrum/fgumi/extract/main'
include { FULCRUM_FGUMI_ALIGN   } from '../modules/fulcrum/fgumi/align/main'
include { FULCRUM_FGUMI_GROUP   } from '../modules/fulcrum/fgumi/group/main'
include { FULCRUM_FGUMI_SIMPLEX } from '../modules/fulcrum/fgumi/simplex/main'
include { FULCRUM_FGUMI_DUPLEX  } from '../modules/fulcrum/fgumi/duplex/main'
include { FULCRUM_FGUMI_FILTER  } from '../modules/fulcrum/fgumi/filter/main'


workflow PREPROCESS_UMI {
    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
    fai                       // channel: [optional] /path/to/reference/fasta_fai, needed for Sentieon
    map_index                 // channel: [mandatory] Pre-computed mapping index
    grouping_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"

    main:

    ch_reports = channel.empty()

    // Extract UMIs from FASTQ, Convert FASTQ files to unmapped BAM with UMI extraction
    FULCRUM_FGUMI_EXTRACT(
        reads
    )

    // Align and sort reads using fgumi fastq + zipper + sort pipeline
    FULCRUM_FGUMI_ALIGN(
        FULCRUM_FGUMI_EXTRACT.out.bam // unmapped bam
    )

    // TODO: Group and call consesnus can also be combined in a single pipe
    // Group reads by UMI
    FULCRUM_FGUMI_GROUP(
        FULCRUM_FGUMI_ALIGN.out.bam // aligned and sorted bam
    )

    // Call consensus reads
    // one process for each mode, gated by meta data is samplesheet?
    FULCRUM_FGUMI_SIMPLEX()
    FULCRUM_FGUMI_DUPLEX()

    ch_consensus_bam = channel.empty()
    ch_consensus_bam = ch_consensus_bam.mix(FULCRUM_FGUMI_SIMPLEX.out.bam)
    ch_consensus_bam = ch_consensus_bam.mix(FULCRUM_FGUMI_DUPLEX.out.bam)

    // TODO: Optional: Duplex metrics
    // FULCRUM_FGUMI_METRICS()

    // Consensus reads are unmapped and must be re-aligned before filtering
    // TODO: Include samtools to allow final outfile to be cram or need to convert out_bam to cram in separate process.
    FULCRUM_FGUMI_FILTER(
        ch_consensus_bam,
        fasta,
        map_index,
    )

    // TODO: Without samtools in FULCRUM_FGUMI_FILTER, we need to convert bam to cram with SAMTOOLS_CONVERT
    // SAMTOOLS_CONVERT(FULCRUM_FGUMI_FILTER.out.bam)

    ch_cram = FULCRUM_FGUMI_FILTER.out.bam

    emit:
    cram = ch_cram // [meta, cram, crai] — duplicate_marked[, recalibrated], indexed
    reports = ch_reports // channel of QC files for MultiQC
}