/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_INTERVALS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prepare genome interval files for scatter/gather operations.

    Takes a BED file of calling regions and produces:
    - A single BED for whole-genome tools
    - A bgzipped + tabixed BED for tools requiring it (Manta)
    - Split BEDs for scatter/gather parallelism (Mutect2)
    - Split bgzipped + tabixed BEDs (Strelka2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_SPLIT } from '../modules/bedtools/split/main'
include { TABIX_BGZIPTABIX as TABIX_ALL } from '../modules/tabix/tabix_bgzip/main'
include { TABIX_BGZIPTABIX as TABIX_SPLIT } from '../modules/tabix/tabix_bgzip/main'

workflow PREPARE_INTERVALS {
    take:
    fai // path: genome .fai index
    intervals_bed // path: BED file of calling regions (or empty)
    _target_bed // path: optional target BED (panel mode, unused for WGS)

    main:


    // Create the intervals channel from the provided BED file
    // If no intervals file is provided, create one from the .fai
    if (intervals_bed) {
        ch_intervals_bed = channel.value(
            [["id": "intervals"], file(intervals_bed, checkIfExists: true)]
        )
    }
    else {
        // Generate a BED from the .fai covering all contigs
        ch_intervals_bed = channel.fromPath(fai, checkIfExists: true)
            .map { fai_file ->
                def bed_content = fai_file
                    .readLines()
                    .collect { line ->
                        def parts = line.split('\t')
                        "${parts[0]}\t0\t${parts[1]}"
                    }
                    .join('\n')
                def bed_file = file("${workDir}/intervals_from_fai.bed")
                bed_file.text = bed_content
                [["id": "intervals"], bed_file]
            }
    }

    // Emit the full intervals BED as-is
    intervals_bed_all = ch_intervals_bed

    // bgzip + tabix the full BED for tools like Manta
    TABIX_ALL(ch_intervals_bed)

    intervals_bed_bgz_tbi_all = TABIX_ALL.out.gz_tbi

    // Split the BED into N chunks for scatter/gather
    BEDTOOLS_SPLIT(
        ch_intervals_bed,
        params.n_interval_splits ?: 10,
    )

    // Flatten the split BEDs into individual [meta, bed] tuples
    intervals_bed_split = BEDTOOLS_SPLIT.out.beds.flatMap { _meta, beds ->
        beds instanceof List
            ? beds.collect { bed -> [["id": bed.baseName], bed] }
            : [[["id": beds.baseName], beds]]
    }

    // bgzip + tabix each split BED for Strelka2
    TABIX_SPLIT(intervals_bed_split)

    intervals_bed_bgz_tbi_split = TABIX_SPLIT.out.gz_tbi

    emit:
    intervals_bed_all = intervals_bed_all
    intervals_bed_bgz_tbi_all = intervals_bed_bgz_tbi_all
    intervals_bed_split = intervals_bed_split
    intervals_bed_bgz_tbi_split = intervals_bed_bgz_tbi_split
}
