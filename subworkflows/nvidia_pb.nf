//
// RUN NVIDIA PARABRICKS
//

include { PARABRICKS_FQ2BAM    } from '../modules/parabricks/fq2bam/main'
include { PARABRICKS_APPLYBQSR } from '../modules/parabricks/applybqsr/main'

workflow PARABRICKS {

    take:
    ch_fastqs
    ch_fasta
    ch_fai
    intervals
    known_sites

    main:

    PARABRICKS_FQ2BAM(
        ch_fastqs,
        ch_fasta,
        ch_fai
    )

    // PARABRICKS_APPLYBQSR

    // CONVERT TO CRAM?

    emit:
    bam = PARABRICKS_FQ2BAM.bam
    bai = PARABRICKS_FQ2BAM.bai
    versions = PARABRICKS_FQ2BAM.versions
}
