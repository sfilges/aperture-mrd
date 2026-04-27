/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALLER_INTERSECTION Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Normalizes SNV VCFs from three somatic callers (Mutect2, Strelka2, LoFreq),
    then intersects using bcftools isec to keep variants with >=2/3 caller agreement.

    SNVs only — indels are excluded from the compendium as they are too noisy for
    MRDetect integration. Indel calling is retained upstream for potential future use.

    Caller-specific handling:
      - Strelka2: SNV VCF used directly (indel VCF excluded)
      - LoFreq:   somatic_final SNV VCF selected (indel VCF excluded)
      - Mutect2:  SNVs extracted via bcftools view --types snps

    Output: Mutect2-annotated VCF of SNVs confirmed by at least one other caller.
    This serves as the patient SNV compendium for downstream MRDetect analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_NORM as NORM_MUTECT2 } from '../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as NORM_STRELKA } from '../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as NORM_LOFREQ } from '../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as NORM_MUSE } from '../modules/bcftools/norm/main'
include { BCFTOOLS_ISEC } from '../modules/bcftools/isec/main'

workflow VCF_CONSENSUS {
    take:
    mutect2_vcf // [meta, filtered.vcf.gz]
    _mutect2_tbi // [meta, filtered.vcf.gz.tbi]
    strelka_snvs_vcf // [meta, somatic_snvs.vcf.gz]
    lofreq_vcf // [meta, [*.vcf.gz]] — multiple LoFreq output VCFs
    muse_vcf // [meta, [*.vcf.gz]] — multiple LoFreq output VCFs
    ch_fasta // [meta, fasta]
    _ch_fasta_fai // [meta, fai]

    main:
    ch_versions = channel.empty()

    //
    // Prepare LoFreq: select somatic_final SNV VCF only (exclude indels and minus-dbsnp)
    //
    ch_lofreq_snvs = lofreq_vcf.map { meta, vcfs ->
        def selected = vcfs instanceof List
            ? vcfs.findAll { f ->
                f.name.endsWith('somatic_final.snvs.vcf.gz')
            }
            : [vcfs]
        [meta, selected.first()]
    }

    //
    // Normalize all three caller VCFs (left-align, split multiallelic)
    //
    NORM_MUTECT2(mutect2_vcf, ch_fasta)
    NORM_STRELKA(strelka_snvs_vcf, ch_fasta)
    NORM_LOFREQ(ch_lofreq_snvs, ch_fasta)
    NORM_MUSE(muse_vcf, ch_fasta)

    ch_versions = ch_versions.mix(NORM_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(NORM_STRELKA.out.versions)
    ch_versions = ch_versions.mix(NORM_LOFREQ.out.versions)
    ch_versions = ch_versions.mix(NORM_MUSE.out.versions)

    //
    // Intersect: keep SNVs present in >=2/3 callers
    // Mutect2 is the first input so -w 1 outputs its annotated records
    //
    ch_isec_input = NORM_MUTECT2.out.vcf
        .join(NORM_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_STRELKA.out.vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_STRELKA.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_LOFREQ.out.vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_LOFREQ.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_MUSE.out.vcf, failOnDuplicate: true, failOnMismatch: true)
        .join(NORM_MUSE.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, m_vcf, m_tbi, s_vcf, s_tbi, l_vcf, l_tbi, mu_vcf, mu_tbi ->
            [meta, [m_vcf, s_vcf, l_vcf, mu_vcf], [m_tbi, s_tbi, l_tbi, mu_tbi]]
        }

    BCFTOOLS_ISEC(ch_isec_input, "+2")
    ch_versions = ch_versions.mix(BCFTOOLS_ISEC.out.versions)

    emit:
    compendium_vcf = BCFTOOLS_ISEC.out.vcf // [meta, compendium.isec.vcf.gz]
    compendium_tbi = BCFTOOLS_ISEC.out.tbi // [meta, compendium.isec.vcf.gz.tbi]
    versions = ch_versions
}
