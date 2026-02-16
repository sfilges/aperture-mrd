


// Mutect2 and post-processing modules
include { GATK4_MUTECT2 as MUTECT2_SOMATIC                              } from '../modules/gatk4/mutect2/main'
include { GATK4_MERGEVCFS as MERGE_MUTECT2                              } from '../modules/gatk4/mergevcfs/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL  } from '../modules/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL         } from '../modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR          } from '../modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION        } from '../modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS as FILTERMUTECTCALLS                  } from '../modules/gatk4/filtermutectcalls/main'

// Other variant callers
include { MANTA_SOMATIC                    } from '../modules/variant_callers/manta/main'
include { STRELKA_SOMATIC                  } from '../modules/variant_callers/strelka2/main'
include { LOFREQ_SOMATIC                   } from '../modules/variant_callers/lofreq/main'
include { MUSE_CALL                        } from '../modules/variant_callers/muse/call/main'

// Utility modules
include { TABIX_BGZIPTABIX                 } from '../modules/tabix/tabix_bgzip/main'
include { BEDTOOLS_SPLIT                   } from '../modules/bedtools/split/main'

workflow TN_SOMATIC_VARIANT_CALLING {

    take:
    cram_variant_calling_pair
    ch_fasta
    ch_fasta_fai
    dict
    germline_resource
    germline_resource_tbi
    versions
    dbsnp
    dbsnp_tbi
    intervals_bed_all
    intervals_bed_gbz_tbi_all
    intervals_bed_split
    intervals_bed_bgz_tbi_split // [meta, intervals.bed, intervals.bed.tbi]


    main:

    //
    // MUTECT2
    //

    // To increase the speed of a Mutect2 workflow, consider running multiple instances 
    // of Mutect2 with e -L (or --intervals) option, with an included BED file of genomic 
    // regions. A target list of 10,000 regions could then be broken up into groups of 100-regions 
    // and run in parallel. In the mutect2 WDL this is achieved by setting the scatter_count parameter. 
    cram_variant_calling_pair_interval_bed_split = cram_variant_calling_pair
        .combine(intervals_bed_split.map { meta, bed -> [ bed ] })

    // Channel: [[meta], normal_cram, normal_crai, tumor_cram, tumor_crai, intervals.bed]
    //cram_variant_calling_pair_interval_bed_split.view()
    
    MUTECT2_SOMATIC(
        cram_variant_calling_pair_interval_bed_split,
        ch_fasta,
        ch_fasta_fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        [],
        []
    )

    // TODO: Add post-processing of Mutect2 output VCF files
    // MERGE VCF files from Mutect2 if multiple intervals were used

    //vcf_to_merge = MUTECT2_SOMATIC.out.vcf
    //MERGE_MUTECT2(vcf_to_merge, dict)

    // GATK best practice is to run the following tools after Mutect2:
    // GATK LearnReadOrientationModel
    //LEARNREADORIENTATIONMODEL()

    //GETPILEUPSUMMARIES_NORMAL()
    //GETPILEUPSUMMARIES_TUMOR()
    // GATK CalculateContamination
    //CALCULATECONTAMINATION()
    // GATK FilterMutectCalls
    //FILTERMUTECTCALLS()
    
    //
    // MANTA and STRELKA
    //

    // Run Manta for somatic SV variant calling and candidate small indels for Strelka
    // Provide all intervals in one file, bg-zipped and tabix indexed
    MANTA_SOMATIC(
        cram_variant_calling_pair,
        ch_fasta,
        ch_fasta_fai,
        intervals_bed_gbz_tbi_all,
        []
    )
    
    // Collect the output candidate small indels VCF and TBI files from Manta for Strelka
    cram_strelka = cram_variant_calling_pair
        .join(MANTA_SOMATIC.out.candidate_small_indels_vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi, failOnDuplicate: true, failOnMismatch: true)
        .combine(intervals_bed_bgz_tbi_split.map { meta, bed, tbi -> [ bed, tbi ] })

    //cram_strelka.view()

    // Strelka calls the entire genome by default, however variant calling may be restricted to an arbitrary 
    // subset of the genome by providing a region file in BED format with the --callRegions configuration option. 
    // The BED file must be bgzip-compressed and tabix-indexed, and only one such BED file may be specified
    STRELKA_SOMATIC(
        cram_strelka,
        ch_fasta,
        ch_fasta_fai
    )

    //
    // LOFREQ
    //

    // TODO: Convert CRAM to BAM for LOFREQ_SOMATIC and MUSE_SOMATIC
    LOFREQ_SOMATIC(
        cram_variant_calling_pair,
        ch_fasta,
        ch_fasta_fai,
        dbsnp,
        dbsnp_tbi
    )

    //
    // MUSE
    //

    //MUSE_SOMATIC(
    //    bam_variant_calling_pair,
    //    ch_fasta
    //)

    //
    // COLLECT OUTPUTS 
    //

    versions = versions.mix(MUTECT2_SOMATIC.out.versions)
    versions = versions.mix(MANTA_SOMATIC.out.versions)
    versions = versions.mix(STRELKA_SOMATIC.out.versions)
    versions = versions.mix(LOFREQ_SOMATIC.out.versions)

    emit:
    versions = versions

}