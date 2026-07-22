//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { SAMTOOLS_FAIDX } from '../modules/samtools/faidx/main'
include { BWA_INDEX } from '../modules/alignment/bwamem/index/main'
include { BWAMEM2_INDEX } from '../modules/alignment/bwamem2/index/main'
include { BWAMEM3_INDEX } from '../modules/alignment/bwamem3/index/main'
include { MINIBWA_INDEX } from '../modules/alignment/minibwa/index/main'

workflow PREPARE_GENOME {
    take:
        fasta_in          // params.fasta
        fai_in            // params.fai
        bwa_index         // params.bwa_index
        bwamem2_index     // params.bwamem2_index
        bwamem3_index     // params.bwamem3_index
        minibwa_index     // params.minibwa_index

    main:

        // Get reference file in fasta format
        ch_fasta = channel.value(
            [["id": "fasta"], file(fasta_in, checkIfExists: true)]
        )
        
        // Generate index if not available
        ch_fai = fai_in
            ? channel.value([["id": "fasta_fai"], file(fai_in, checkIfExists: true)])
            : SAMTOOLS_FAIDX(ch_fasta).fai.first()

        // Generate index for chosen aligner if not available
        if (params.aligner == 'bwamem3') {
            ch_index = bwamem3_index
                ? channel.value([[id:'bwamem3_index'], file("${bwamem3_index}/*", checkIfExists: true)])
                : BWAMEM3_INDEX(ch_fasta).index.first()
        } else if (params.aligner == 'minibwa') {
            ch_index = minibwa_index
                ? channel.value([[id:'minibwa_index'], file("${minibwa_index}/*", checkIfExists: true)])
                : MINIBWA_INDEX(ch_fasta).index.first()
        } else if (params.aligner == 'bwamem') {
            ch_index = bwa_index
                ? channel.value([[id:'bwa_index'], file("${bwa_index}/*", checkIfExists: true)])
                : BWA_INDEX(ch_fasta).index.first()
        } else if (params.aligner == 'bwamem2') {
            ch_index = bwamem2_index
                ? channel.value([[id:'bwamem2_index'], file("${bwamem2_index}/*", checkIfExists: true)])
                : BWAMEM2_INDEX(ch_fasta).index.first()
        } else {
            error "PREPARE_GENOME: unsupported aligner '${params.aligner}'"
        }
    
    emit:
        fasta = ch_fasta
        fai   = ch_fai
        index = ch_index
}
