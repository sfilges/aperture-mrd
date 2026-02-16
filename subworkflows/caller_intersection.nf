/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALLER_INTERSECTION Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Normalizes VCFs from three somatic callers (Mutect2, Strelka2, LoFreq),
    then intersects using bcftools isec to keep variants with >=2/3 caller agreement.

    Caller-specific handling:
      - Strelka2: SNV and indel VCFs are concatenated into a single VCF
      - LoFreq:   somatic_final SNV + indel VCFs are selected and concatenated
      - Mutect2:  used as-is (already a single filtered VCF)

    Output: Mutect2-annotated VCF of variants confirmed by at least one other caller.
    This serves as the patient SNV compendium for downstream MRDetect analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_CONCAT as CONCAT_STRELKA } from '../modules/bcftools/concat/main'
include { BCFTOOLS_CONCAT as CONCAT_LOFREQ  } from '../modules/bcftools/concat/main'
include { BCFTOOLS_NORM as NORM_MUTECT2     } from '../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as NORM_STRELKA     } from '../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as NORM_LOFREQ      } from '../modules/bcftools/norm/main'
include { BCFTOOLS_ISEC                     } from '../modules/bcftools/isec/main'

workflow CALLER_INTERSECTION {

    take:
    mutect2_vcf        // [meta, filtered.vcf.gz]
    mutect2_tbi        // [meta, filtered.vcf.gz.tbi]
    strelka_snvs_vcf   // [meta, somatic_snvs.vcf.gz]
    strelka_indels_vcf // [meta, somatic_indels.vcf.gz]
    lofreq_vcf         // [meta, [*.vcf.gz]] — multiple LoFreq output VCFs
    ch_fasta           // [meta, fasta]
    ch_fasta_fai       // [meta, fai]

    main:
    ch_versions = channel.empty()

    //
    // Prepare Strelka2: concatenate SNV + indel VCFs into a single VCF
    //
    ch_strelka_to_concat = strelka_snvs_vcf
        .join(strelka_indels_vcf, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, snvs, indels -> [meta, [snvs, indels]] }

    CONCAT_STRELKA(ch_strelka_to_concat)
    ch_versions = ch_versions.mix(CONCAT_STRELKA.out.versions)

    //
    // Prepare LoFreq: select somatic_final SNV + indel VCFs (exclude minus-dbsnp duplicates)
    //
    ch_lofreq_to_concat = lofreq_vcf.map { meta, vcfs ->
        def selected = vcfs instanceof List
            ? vcfs.findAll { f ->
                f.name.endsWith('somatic_final.snvs.vcf.gz') ||
                f.name.endsWith('somatic_final.indels.vcf.gz')
            }
            : [vcfs]
        [meta, selected]
    }

    CONCAT_LOFREQ(ch_lofreq_to_concat)
    ch_versions = ch_versions.mix(CONCAT_LOFREQ.out.versions)

    //
    // Normalize all three caller VCFs (left-align, split multiallelic)
    //
    NORM_MUTECT2(mutect2_vcf, ch_fasta)
    NORM_STRELKA(CONCAT_STRELKA.out.vcf, ch_fasta)
    NORM_LOFREQ(CONCAT_LOFREQ.out.vcf, ch_fasta)

    ch_versions = ch_versions.mix(NORM_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(NORM_STRELKA.out.versions)
    ch_versions = ch_versions.mix(NORM_LOFREQ.out.versions)

    //
    // Intersect: keep variants present in >=2/3 callers
    // Mutect2 is the first input so -w 1 outputs its annotated records
    //
    ch_isec_input = NORM_MUTECT2.out.vcf
        .join(NORM_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_STRELKA.out.vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_STRELKA.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_LOFREQ.out.vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_LOFREQ.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, m_vcf, m_tbi, s_vcf, s_tbi, l_vcf, l_tbi ->
            [meta, [m_vcf, s_vcf, l_vcf], [m_tbi, s_tbi, l_tbi]]
        }

    BCFTOOLS_ISEC(ch_isec_input, "+2")
    ch_versions = ch_versions.mix(BCFTOOLS_ISEC.out.versions)

    emit:
    compendium_vcf = BCFTOOLS_ISEC.out.vcf  // [meta, compendium.isec.vcf.gz]
    compendium_tbi = BCFTOOLS_ISEC.out.tbi  // [meta, compendium.isec.vcf.gz.tbi]
    versions       = ch_versions
}
