//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { SAMTOOLS_FAIDX } from '../modules/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
        reference

    main:
        // checks the reference and its indexes, if the indexes are not there creates them
        reference_file = file(reference)
        if (reference_file.isEmpty()) {
            log.error "--reference points to a non existing file"
            exit 1
        }


    emit:
        checked_reference = reference
}