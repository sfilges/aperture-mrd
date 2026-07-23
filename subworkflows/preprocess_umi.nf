/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMI Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fulcrum fgumi implementation (Rust-based, supersedes previous fgbio)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FULCRUM_FGUMI_EXTRACT } from '../modules/fulcrum/fgumi/extract/main'
include { FULCRUM_FGUMI_ALIGN   } from '../modules/fulcrum/fgumi/align/main'
include { FULCRUM_FGUMI_GROUPREADSBYUMI } from '../modules/fulcrum/fgumi/groupreadsbyumi/main'
include { FULCRUM_FGUMI_SIMPLEX } from '../modules/fulcrum/fgumi/simplex/main'
include { FULCRUM_FGUMI_DUPLEX } from '../modules/fulcrum/fgumi/duplex/main'
include { FULCRUM_FGUMI_FILTER } from '../modules/fulcrum/fgumi/filter/main'


workflow FASTQ_CREATE_UMI_CONSENSUS_FGBIO {
    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
    fai                       // channel: [optional] /path/to/reference/fasta_fai, needed for Sentieon
    map_index                 // channel: [mandatory] Pre-computed mapping index
    groupreadsbyumi_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"

    main:

    // Extract UMIs from FASTQ:
    FULCRUM_FGUMI_EXTRACT(reads)

    // Align and sort reads using fgumi fastq + zipper + sort pipeline
    FULCRUM_FGUMI_ALIGN(
        FULCRUM_FGUMI_EXTRACT.out.bam
    )

    // Group reads by UMI
    FULCRUM_FGUMI_GROUPREADSBYUMI()

    // Call consensus reads
    // one process for each mode, gated by meta data is samplesheet?
    FULCRUM_FGUMI_SIMPLEX()
    FULCRUM_FGUMI_DUPLEX()

    // Optional: Duplex metrics

    // Filter
    FULCRUM_FGUMI_FILTER()

}