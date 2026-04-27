#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Stefan Filges / Aperture-MRD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Started July 2024.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sfilges/aperture-mrd:
        An  analysis pipeline to detect somatic variants for tumor-informed
        circulating tumor DNA detection from whole genome sequencing.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sfilges/aperture-mrd
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Enable dsl 2
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from './subworkflows/utils'
include { PREPARE_INTERVALS } from './subworkflows/prepare_intervals'
include { PREPROCESS_GATK } from './subworkflows/preprocess_gatk'
include { TN_SOMATIC_VARIANT_CALLING } from './subworkflows/tn_somatic_variant_calling'
include { VCF_CONSENSUS } from './subworkflows/vcf_consensus'
include { VCF_FILTER } from './subworkflows/vcf_filter'
include { MULTIQC } from './modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {

    if (params.validate_params) {
        // Validate parameters relative to the parameter JSON schema 
        // in default location: "nextflow_schema.json"
        validateParameters()
    }

    // Initilialise versions and reports channels
    ch_versions = channel.empty()
    _ch_reports = channel.empty()

    ch_multiqc_files = channel.empty()

    // It is the order of fields in the samplesheet JSON schema which defines 
    // the order of items in the channel, *not* the order of fields in the 
    // samplesheet file.
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .map { meta, fastq1, fastq2 ->
            // structure the output depending on the input
            if (fastq2) {
                [meta, [fastq1, fastq2]]
            }
            else if (fastq1) {
                [meta, [fastq1]]
            }
        }

    //ch_samplesheet.view()

    //
    // Prepare all input files
    //

    // TODO Add prepare genomes subworkflow to prepare reference files 
    // to have them available for all downstream processes.

    // Get reference file in fasta format
    ch_fasta = channel.value(
        [["id": "fasta"], file(params.fasta, checkIfExists: true)]
    )

    ch_fasta_fai = channel.value(
        [["id": "fasta_fai"], file("${params.fasta}.fai", checkIfExists: true)]
    )

    // Create known sites channels (collect once, reuse multiple times)
    known_indels = params.known_indels ? channel.fromPath(params.known_indels).collect() : channel.value([])
    known_indels_tbi = params.known_indels_tbi ? channel.fromPath(params.known_indels_tbi).collect() : channel.value([])
    dbsnp = params.dbsnp ? channel.fromPath(params.dbsnp).collect() : channel.value([])
    dbsnp_tbi = params.dbsnp ? channel.fromPath("${params.dbsnp}.tbi").collect() : channel.value([])
    _known_snps = params.known_snps ? channel.fromPath(params.known_snps).collect() : channel.value([])
    _known_snps_tbi = params.known_snps ? channel.fromPath("${params.known_snps}.tbi").collect() : channel.value([])
    germline_resource = params.germline_resource ? channel.fromPath(params.germline_resource).collect() : channel.value([])
    germline_resource_tbi = params.germline_resource_tbi ? channel.fromPath("${params.germline_resource_tbi}").collect() : channel.value([])

    // TODO Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals

    PREPARE_INTERVALS(
        params.fai,
        params.intervals,
        [],
    )

    // Combine known sites and dict — needed by PREPROCESS_READS and variant calling
    known_sites_indels = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()
    dict = channel.fromPath(params.dict, checkIfExists: true).collect()

    //
    // Preprocessing: FastQC → fastp → alignment → MarkDuplicates → BQSR → QC
    //
    PREPROCESS_GATK(
        ch_samplesheet,
        ch_fasta,
        ch_fasta_fai,
        dict,
        known_sites_indels,
        known_sites_indels_tbi,
    )

    ch_cram_for_variant_calling = PREPROCESS_GATK.out.cram
    ch_versions = ch_versions.mix(PREPROCESS_GATK.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_GATK.out.reports)

    //
    // Logic to combine tumor-normal pairs. Does *not* work for tumor only or germline only samples!
    //

    //The branch operator forwards each item from a source channel to one of multiple 
    //output channels, based on a selection criteria.
    ch_cram_variant_calling_by_status = ch_cram_for_variant_calling.branch { row ->
        normal: row[0].status == 0
        tumor: row[0].status == 1
    }

    // All Germline samples
    cram_variant_calling_normal_to_cross = ch_cram_variant_calling_by_status.normal.map { meta, cram, crai -> [meta.sample, meta, cram, crai] }

    // All tumor samples
    cram_variant_calling_pair_to_cross = ch_cram_variant_calling_by_status.tumor.map { meta, cram, crai -> [meta.sample, meta, cram, crai] }

    // Tumor - normal pairs
    // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
    cram_variant_calling_pair = cram_variant_calling_normal_to_cross
        .cross(cram_variant_calling_pair_to_cross)
        .map { normal, tumor ->
            def meta = [:]

            meta.id = "${tumor[1].id}_vs_${normal[1].id}".toString()
            meta.normal_id = normal[1].id
            meta.sample = normal[0]
            meta.tumor_id = tumor[1].id

            [meta, normal[2], normal[3], tumor[2], tumor[3]]
        }

    // [meta, normal_cram, normal_crai, tumor_cram, tumor_crai]
    //cram_variant_calling_pair.view()

    //
    // Run SNV variant calling on tumor-normal pairs
    //

    TN_SOMATIC_VARIANT_CALLING(
        cram_variant_calling_pair,
        ch_fasta,
        ch_fasta_fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        ch_versions,
        dbsnp,
        dbsnp_tbi,
        PREPARE_INTERVALS.out.intervals_bed_all,
        PREPARE_INTERVALS.out.intervals_bed_bgz_tbi_all,
        PREPARE_INTERVALS.out.intervals_bed_split,
        PREPARE_INTERVALS.out.intervals_bed_bgz_tbi_split,
    )

    //
    // Caller intersection: keep SNVs with >=2/3 caller agreement
    //
    VCF_CONSENSUS(
        TN_SOMATIC_VARIANT_CALLING.out.mutect2_vcf,
        TN_SOMATIC_VARIANT_CALLING.out.mutect2_tbi,
        TN_SOMATIC_VARIANT_CALLING.out.strelka_snvs_vcf,
        TN_SOMATIC_VARIANT_CALLING.out.lofreq_vcf,
        ch_fasta,
        ch_fasta_fai,
    )

    ch_versions = ch_versions.mix(VCF_CONSENSUS.out.versions)

    //
    // Blacklist filtering + gnomAD common variant exclusion
    //
    ch_blacklists = channel.fromPath(
            [
                params.encode_blacklist,
                params.centromeres,
                params.simple_repeats,
            ]
        )
        .collect()

    VCF_FILTER(
        VCF_CONSENSUS.out.compendium_vcf,
        VCF_CONSENSUS.out.compendium_tbi,
        germline_resource,
        germline_resource_tbi,
        ch_blacklists,
    )
    ch_versions = ch_versions.mix(VCF_FILTER.out.versions)

    // TODO: CNA calling with CNVkit



    // TODO: Annotate variants with VEP


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/${workflow.runName}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
        .set { _ch_collated_versions }

    //
    // Multiqc
    //
    MULTIQC(
        ch_multiqc_files.collect()
    )
}
