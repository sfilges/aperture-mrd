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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters              } from 'plugin/nf-schema'
include { samplesheetToList               } from 'plugin/nf-schema'
include { PREPARE_GENOME                  } from './subworkflows/prepare_genome'
include { PREPARE_INTERVALS               } from './subworkflows/prepare_intervals'
include { PREPROCESS_READS                } from './subworkflows/preprocess_reads'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from './modules/samtools/convert/main'
include { TN_SOMATIC_SNV_CALLING          } from './subworkflows/tn_somatic_snv_calling'
include { TN_SOMATIC_SIGNATURES           } from './subworkflows/tn_somatic_signatures'
include { VCF_CONSENSUS                   } from './subworkflows/vcf_consensus'
include { VCF_FILTER                      } from './subworkflows/vcf_filter'
include { MULTIQC                         } from './modules/multiqc/main'
include { pairTumorNormal                 } from './subworkflows/utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    if (params.validate_params) {
        // Validate parameters relative to the parameter JSON schema in default location: "nextflow_schema.json"
        // As of 2026-07-23 no automatic builder for samplesheet schema exist, you must build your own from the example.
        // See samplesheet schema example: https://nextflow-io.github.io/nf-schema/latest/nextflow_schema/sample_sheet_schema_examples/
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

    // Panel of normals for variant calling (uses 1000 Genomes by default)
    pon = params.pon ? channel.fromPath("${params.pon}").collect() : channel.value([])
    pon_tbi = params.pon_tbi ? channel.fromPath("${params.pon_tbi}").collect() : channel.value([])

    // Decide which intervals to use. PREPARE_INTERVALS resolves this itself with
    // file(), so it must be a path (or null) — passing a channel here makes that
    // call throw. A null falls back to building intervals from the .fai.
    def intervals = null
    if (mode == 'wes') {
        intervals = params.wes_use_baits_as_intervals ? params.wes_baits : params.wes_targets
    }
    else if (mode == 'wgs') {
        intervals = params.wgs_intervals
    }
    else {
        error("Unknown mode: ${mode} (expected 'wes' or 'wgs')")
    }

    // TODO Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals
    PREPARE_INTERVALS(
        PREPARE_GENOME.out.fai,
        intervals,
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

    // TODO: If preprocessing == 'umi', use the UMI workflow instead
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
    // CRAM to BAM conversion for tools that cannot (efficiently) read CRAM:
    // CNVkit, MSIsensor2, MuSE, LoFreq. Conversion happens per sample, before
    // pairing, so a normal shared by several tumors is only converted once
    // (same as: https://github.com/nf-core/sarek/blob/3a6a502a93c3e19d3699fa1be682e801bf2745ad/workflows/sarek/main.nf#L33).
    // The BAMs are unpublished intermediates; CRAM stays the archival format.
    //
    CRAM_TO_BAM(
        ch_cram_for_variant_calling.map { meta, cram, _crai -> [meta, cram] },
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
    )
    ch_bam_for_variant_calling = CRAM_TO_BAM.out.bam
        .join(CRAM_TO_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)

    //
    // Combine tumor-normal pairs, in parallel for CRAM and BAM.
    // Both channels: [meta, normal_file, normal_index, tumor_file, tumor_index]
    //
    cram_variant_calling_pair = pairTumorNormal(ch_cram_for_variant_calling)
    bam_variant_calling_pair = pairTumorNormal(ch_bam_for_variant_calling)

    //
    // Run SNV variant calling on tumor-normal pairs
    //

    TN_SOMATIC_SNV_CALLING(
        cram_variant_calling_pair,
        bam_variant_calling_pair,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        dbsnp,
        dbsnp_tbi,
        pon,
        pon_tbi,
        PREPARE_INTERVALS.out.intervals_bed_all,
        PREPARE_INTERVALS.out.intervals_bed_bgz_tbi_all,
        PREPARE_INTERVALS.out.intervals_bed_split,
        PREPARE_INTERVALS.out.intervals_bed_bgz_tbi_split,
    )

    //
    // Caller intersection: keep SNVs with >=2 caller agreement
    //
    VCF_CONSENSUS(
        TN_SOMATIC_SNV_CALLING.out.mutect2_vcf,
        TN_SOMATIC_SNV_CALLING.out.mutect2_tbi,
        TN_SOMATIC_SNV_CALLING.out.strelka_snvs_vcf,
        //TN_SOMATIC_SNV_CALLING.out.muse_vcf,
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

    // TODO: CNA calling with CNVkit — bam_variant_calling_pair matches the
    // SOMATIC_CNV_CALLING take order [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]

    // TODO: Mutational signature detection. MSIsensor2 is tumor-only, so feed it
    // from the per-sample BAM channel to keep the original sample meta:
    TN_SOMATIC_SIGNATURES(
        ch_bam_for_variant_calling.filter { meta, _bam, _bai -> meta.status == 1 },
        params.msisensor2_models
    )

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
        .collectFile(storeDir: "${params.outdir}/${params.run_id}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
        .set { _ch_collated_versions }

    //
    // Multiqc
    //
    MULTIQC(
        ch_multiqc_files.collect()
    )
}
