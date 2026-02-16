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

include { PREPARE_INTERVALS } from './subworkflows/prepare_intervals/main'

// Load preprocessing modules
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc/main'
include { FASTP } from './modules/fastp/main'
include { BWA_INDEX } from './modules/bwamem/index/main'
include { BWA_MEM } from './modules/bwamem/mem/main'
include { GATK4_MARKDUPLICATES } from './modules/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from './modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from './modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from './modules/samtools/index/main'
//include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from './modules/samtools/convert/main'
include { BWAMEM2_INDEX } from './modules/bwamem2/index/main'
include { BWAMEM2_MEM } from './modules/bwamem2/mem/main'


// Load parabricks modules
//include { PARABRICKS_FQ2BAM                 } from './modules/parabricks/fq2bam/main'
//include { PARABRICKS_APPLYBQSR              } from './modules/parabricks/applybqsr/main'
//include { PARABRICKS_COLLECTMULTIPLEMETRICS } from './modules/parabricks/collectmultiplemetrics/main'
//include { PARABRICKS_MUTECTCALLER           } from './modules/parabricks/mutectcaller/main'

// Generate report metrics
//include { GATK4_COLLECTMULTIPLEMETRICS } from './modules/gatk4/collectmultiplemetrics/main'
// include { GATK4_COLLECTWGSMETRICS    } from './modules/gatk4/collectwgsmetrics/main'
// include { NGSCHECKMATE               } from './modules/ngscheckmate/main'
include { SAMTOOLS_STATS } from './modules/samtools/stats/main'
include { MOSDEPTH } from './modules/mosdepth/main'
include { MULTIQC } from './modules/multiqc/main'

include { TN_SOMATIC_VARIANT_CALLING } from './subworkflows/tn_somatic_variant_calling'

// Variant annotation


// Variant merging & filtering


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
    ch_versions = Channel.empty()
    ch_reports = Channel.empty()

    ch_multiqc_files = Channel.empty()

    // It is the order of fields in the samplesheet JSON schema which defines 
    // the order of items in the channel, *not* the order of fields in the 
    // samplesheet file.
    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
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
    ch_fasta = Channel.value(
        [["id": "fasta"], file(params.fasta, checkIfExists: true)]
    )

    ch_fasta_fai = Channel.value(
        [["id": "fasta_fai"], file("${params.fasta}.fai", checkIfExists: true)]
    )

    // Create known sites channels (collect once, reuse multiple times)
    known_indels = params.known_indels ? Channel.fromPath(params.known_indels).collect() : Channel.value([])
    known_indels_tbi = params.known_indels_tbi ? Channel.fromPath(params.known_indels_tbi).collect() : Channel.value([])
    dbsnp = params.dbsnp ? Channel.fromPath(params.dbsnp).collect() : Channel.value([])
    dbsnp_tbi = params.dbsnp ? Channel.fromPath("${params.dbsnp}.tbi").collect() : Channel.value([])
    known_snps = params.known_snps ? Channel.fromPath(params.known_snps).collect() : Channel.value([])
    known_snps_tbi = params.known_snps ? Channel.fromPath("${params.known_snps}.tbi").collect() : Channel.value([])
    germline_resource = params.germline_resource ? Channel.fromPath(params.germline_resource).collect() : Channel.value([])
    germline_resource_tbi = params.germline_resource_tbi ? Channel.fromPath("${params.germline_resource_tbi}").collect() : Channel.value([])

    // TODO Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals

    PREPARE_INTERVALS(
        params.fai,
        params.intervals,
        [],
    )

    //
    // Preprocessing
    //

    // TODO Merge into a preprocessing subworkflow
    FASTQC_RAW(ch_samplesheet)
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

    ch_fastqc_raw_html = FASTQC_RAW.out.html
    ch_fastqc_raw_zip = FASTQC_RAW.out.zip

    FASTP(
        ch_samplesheet,
        params.merge_fastqs,
        params.split_fastq,
        params.trim_nextseq,
        params.length_required,
        params.save_fastqs,
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    if (params.split_fastq) {
        ch_trim_reads = FASTP.out.reads
            .map { meta, reads ->
                def read_files = reads.sort(false) { a, b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                [meta + [n_fastq: read_files.size()], read_files]
            }
            .transpose()
    }
    else {
        ch_trim_reads = FASTP.out.reads
    }

    if (params.merge_fastqs) {
        ch_trim_reads = ch_trim_reads.mix(FASTP.out.reads_merged)
    }

    //FASTQC_TRIM(ch_trim_reads)

    //ch_trim_reads.view()

    // Output file channels and collect for multiqc
    ch_trim_json = FASTP.out.json
    //ch_fastqc_trim_html = FASTQC_TRIM.out.html
    //ch_fastqc_trim_zip  = FASTQC_TRIM.out.zip

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_fastqc_raw_zip.collect { it[1] }.ifEmpty([])
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_trim_json.collect { it[1] }.ifEmpty([])
    )

    //ch_multiqc_files = ch_multiqc_files.mix(
    //    ch_fastqc_trim_zip.collect{it[1]}.ifEmpty([])
    //)


    //
    // Mapping
    //

    // ALIGNMENT
    // Currently hardcoded to add read groups 
    // https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
    // Based on official SAM specification: https://samtools.github.io/hts-specs/


    if (params.aligner == 'bwa') {
        if (params.bwaindex) {
            // Create a channel from existing index files
            ch_index = Channel.value(
                [["id": "bwa_index"], file("${params.bwaindex}/*{,.amb,.ann,.bwt,.pac,.sa}")]
            )
        }
        else {
            // Generate index with BWA_INDEX process
            BWA_INDEX(ch_fasta)
            ch_index = BWA_INDEX.out.index
        }

        BWA_MEM(
            ch_trim_reads,
            ch_index,
        )

        ch_aligned_bams = BWA_MEM.out.bam
    }
    else if (params.aligner == 'bwamem2') {
        if (params.bwamem2index) {
            // Create a channel from existing index files
            ch_index = Channel.value(
                [["id": "bwamem2_index"], file("${params.bwamem2index}/*{,.amb,.ann,.bwt,.pac,.sa}")]
            )
        }
        else {
            // Generate index with BWAMEM2_INDEX process
            BWAMEM2_INDEX(ch_fasta)
            ch_index = BWAMEM2_INDEX.out.index
        }

        BWAMEM2_MEM(
            ch_trim_reads,
            ch_fasta,
            ch_index,
        )

        ch_aligned_bams = BWAMEM2_MEM.out.bam
    }
    else {
        error("Unknown aligner: ${params.aligner}. Supported aligners are: bwa, bwamem2.")
    }

    //ch_index.view()

    // By default, if you don’t specify a size, the groupTuple operator will not 
    // emit any groups until all inputs have been received. If possible, you should 
    // always try to specify the number of expected elements in each group using the 
    // size option, so that each group can be emitted as soon as it’s ready. 
    // In cases where the size of each group varies based on the grouping key, 
    // you can use the built-in groupKey() function, which allows you to define a 
    // different expected size for each group.
    bam_grouped = ch_aligned_bams
        .map { meta, bam ->
            [groupKey(meta, meta.n_fastq), bam]
        }
        .groupTuple()
    // the first element of each tuple is used as the grouping key.

    //bam_grouped.view()

    GATK4_MARKDUPLICATES(
        bam_grouped,
        ch_fasta,
        ch_fasta_fai,
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    // Combine known sites
    known_sites_indels = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    dict = Channel.fromPath(params.dict, checkIfExists: true).collect()

    // Ensure proper channel broadcasting for BASERECALIBRATOR
    GATK4_BASERECALIBRATOR(
        GATK4_MARKDUPLICATES.out.cram,
        ch_fasta,
        ch_fasta_fai,
        dict,
        known_sites_indels,
        known_sites_indels_tbi,
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    ch_multiqc_files = ch_multiqc_files.mix(GATK4_BASERECALIBRATOR.out.table.collect { meta, table -> [table] })

    // Collect the output table for BQSR
    cram_applybqsr = GATK4_MARKDUPLICATES.out.cram.join(GATK4_BASERECALIBRATOR.out.table, failOnDuplicate: true, failOnMismatch: true)
    //cram_applybqsr.view()

    // Apply BQSR
    GATK4_APPLYBQSR(
        cram_applybqsr,
        ch_fasta,
        ch_fasta_fai,
        dict,
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    INDEX_CRAM(
        GATK4_APPLYBQSR.out.cram
    )

    // Join with the crai file, change output tuple to [meta, cram, crai]
    ch_cram_for_variant_calling = GATK4_APPLYBQSR.out.cram
        .join(
            INDEX_CRAM.out.crai,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { meta, cram, crai -> [meta, cram, crai] }

    // Collect Alignment metrics

    // CollectMultipleMetrics combines multiple GATK4 metrics collection tools into one process.
    // This reduces the number of processes and allows for more efficient resource usage.
    // The downside is that it is currently not be as flexible as using individual tools.
    //GATK4_COLLECTMULTIPLEMETRICS(GATK4_APPLYBQSR.out.cram, ch_fasta, ch_fasta_fai)

    SAMTOOLS_STATS(ch_cram_for_variant_calling, ch_fasta)

    MOSDEPTH(ch_cram_for_variant_calling, ch_fasta)

    // Collect the output metrics files for MultiQC
    //[
    //    GATK4_COLLECTMULTIPLEMETRICS.out.alignment_summary_metrics,
    //    GATK4_COLLECTMULTIPLEMETRICS.out.insert_size_metrics,
    //    GATK4_COLLECTMULTIPLEMETRICS.out.quality_by_cycle_metrics,
    //    GATK4_COLLECTMULTIPLEMETRICS.out.base_distribution_by_cycle_metrics,
    //    GATK4_COLLECTMULTIPLEMETRICS.out.gc_bias_metrics,
    //    GATK4_COLLECTMULTIPLEMETRICS.out.quality_yield_metrics
    //].each { ch ->
    //    ch_multiqc_files = ch_multiqc_files.mix(ch)
    //}

    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { meta, stats -> [stats] })
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect { meta, global_txt -> [global_txt] })
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.collect { meta, regions_txt -> [regions_txt] })

    // Gather versions of all tools used
    //ch_versions = ch_versions.mix(GATK4_COLLECTMULTIPLEMETRICS.out.versions)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    // BAM_NGSCHECKMATE



    //
    // Logic to combine tumor-normal pairs. Does *not* work for tumor only or germline only samples!
    //

    //The branch operator forwards each item from a source channel to one of multiple 
    //output channels, based on a selection criteria.
    ch_cram_variant_calling_by_status = ch_cram_for_variant_calling.branch {
        normal: it[0].status == 0
        tumor: it[0].status == 1
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

    // VARIANT QC

    // VARIANT ANNOTATION

    // VARIANT MERGING & FILTERING

    // 1. Intersect variants common betwen at least 2 callers

    // 2. Apply blacklist filtering on intersected variants
    // Maybe the lists should be merged into a single blacklist file?
    // 2.1 ENCODE blacklist
    // 2.2 gnomad blacklist
    // 2.3 Repeat regions
    // 2.4 Other blacklists



    // 3. Apply additional filters



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/${workflow.runName}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // Multiqc
    //
    MULTIQC(
        ch_multiqc_files.collect()
    )
}
