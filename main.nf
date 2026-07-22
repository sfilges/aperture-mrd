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
include { PREPARE_GENOME } from './subworkflows/prepare_genome'
include { PREPARE_INTERVALS } from './subworkflows/prepare_intervals'
include { PREPROCESS_READS } from './subworkflows/preprocess_reads'
include { TN_SOMATIC_SNV_CALLING } from './subworkflows/tn_somatic_snv_calling'
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

    // Initilialise reports channel
    _ch_reports = channel.empty()

    ch_multiqc_files = channel.empty()

    // It is the order of fields in the samplesheet JSON schema which defines 
    // the order of items in the channel, *not* the order of fields in the 
    // samplesheet file.
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .map { meta, fastq1, fastq2 ->
            def CN = params.seq_center ? "CN:${params.seq_center}\t" : ''
            meta = meta + [read_group: "@RG\tID:${meta.id}\t${CN}SM:${meta.id}\tLB:${meta.id}\tPL:${params.seq_platform}"]
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
    PREPARE_GENOME(
        params.fasta,
        params.fai,
        params.bwa_index,
        params.bwamem2_index,
        params.bwamem3_index,
        params.minibwa_index,
    )


    // Create WES file channels
    mode = params.mode
    wes_baits = params.wes_baits ? channel.fromPath(params.wes_baits, checkIfExists: true).collect() : channel.value([])
    wes_targets = params.wes_targets ? channel.fromPath(params.wes_targets, checkIfExists: true).collect() : channel.value([])

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
        PREPARE_GENOME.out.fai,
        params.intervals,
        [],
    )

    // Combine known sites and dict — needed by PREPROCESS_READS (gatk mode) and variant calling
    known_sites_indels = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()
    dict = channel.fromPath(params.dict, checkIfExists: true).collect()

    //
    // Preprocessing: fast or gatk
    //
    if (!(params.preprocessing in ['fast', 'gatk'])) {
        error("Unknown preprocessing mode: ${params.preprocessing}")
    }

    if (params.preprocessing == "fast" && params.aligner in ['bwamem', 'bwamem2']) {
        log.warn "preprocessing='fast' with aligner='${params.aligner}': fast mode is tuned for bwamem3/minibwa; '${params.aligner}' will run but may be suboptimal."
    }

    // FastQC → fastp → alignment → MarkDuplicates → [BQSR] → indexing → QC
    PREPROCESS_READS(
        ch_samplesheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.index,
        mode,
        wes_baits,
        wes_targets,
        dict,
        known_sites_indels,
        known_sites_indels_tbi,
    )

    ch_cram_for_variant_calling = PREPROCESS_READS.out.cram
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_READS.out.reports)

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

    TN_SOMATIC_SNV_CALLING(
        cram_variant_calling_pair,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        dict,
        germline_resource,
        germline_resource_tbi,
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
        TN_SOMATIC_SNV_CALLING.out.mutect2_vcf,
        TN_SOMATIC_SNV_CALLING.out.mutect2_tbi,
        TN_SOMATIC_SNV_CALLING.out.strelka_snvs_vcf,
        TN_SOMATIC_SNV_CALLING.out.lofreq_vcf,
        TN_SOMATIC_SNV_CALLING.out.muse_vcf,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
    )

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

    // TODO: CNA calling with CNVkit



    // TODO: Annotate variants with VEP


    //
    // Collate and save software versions
    //
    // Versions are emitted by each process onto the global "versions" topic
    // channel as (process, tool, version) tuples (Nextflow >=25.04).
    channel.topic('versions')
        .distinct()
        .map { process, tool, version ->
            [process.tokenize(':').last(), "  ${tool}: ${version}"]
        }
        .groupTuple()
        .map { process, tool_versions ->
            "\"${process}\":\n${tool_versions.unique().sort().join('\n')}"
        }
        .collectFile(storeDir: "${params.outdir}/${workflow.runName}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
        .set { _ch_collated_versions }

    //
    // Multiqc
    //
    MULTIQC(
        ch_multiqc_files.collect()
    )
}
