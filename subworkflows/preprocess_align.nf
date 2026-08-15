/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGNMENT Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs only one aligner set by params.aligner using ext.when as defined in
    conf/modules.config.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { BWA_MEM     } from '../modules/alignment/bwamem/mem/main'
include { BWAMEM2_MEM } from '../modules/alignment/bwamem2/mem/main'
include { BWAMEM3_MEM } from '../modules/alignment/bwamem3/mem/main'
include { MINIBWA_MAP } from '../modules/alignment/minibwa/map/main'


workflow PREPROCESS_ALIGN {
    take:
        ch_trim_reads
        ch_fasta
        ch_fai
        ch_index

    main:
    // Only one aligner is run based on ext.when
    BWA_MEM(
        ch_trim_reads,
        ch_index,
        ch_fasta,
    )

    BWAMEM2_MEM(
        ch_trim_reads,
        ch_fasta,
        ch_index,
    )

    BWAMEM3_MEM(
        ch_trim_reads,
        ch_index,
        ch_fasta,
    )

    MINIBWA_MAP(
        ch_trim_reads,
        ch_index,
        ch_fasta,
        ch_fai,
    )

    // Get the bam files from the aligner, only one has values
    ch_aligned_crams = channel.empty()
    ch_aligned_crams = ch_aligned_crams.mix(BWA_MEM.out.cram)
    ch_aligned_crams = ch_aligned_crams.mix(BWAMEM2_MEM.out.cram)
    ch_aligned_crams = ch_aligned_crams.mix(BWAMEM3_MEM.out.cram)
    ch_aligned_crams = ch_aligned_crams.mix(MINIBWA_MAP.out.cram)

    emit:
        crams = ch_aligned_crams
}