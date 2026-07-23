/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT ANNOTATION Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Annotates VCF files using Enembl VEP.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ENSEMBLVEP_VEP } from '../modules/annotation/ensemblvep/vep/main'
include { FASTVEP_ANNOTATE } from '../modules/annotation/fastvep/annotate/main' 

workflow VCF_VEP_ANNOTATE {
    take:
        vcf_for_vep           // tuple val(meta), path(vcf)
        vep_genome        // params.vep_genome
        vep_species       // params.vep_species
        vep_cache_version // params.vep_cache_version
        vep_cache         // params.vep_cache
        fasta         // tuple val(meta2), path(fasta)

    main:
    // TODO: Need to download cache first if not available (in separate workflow?)
    // TODO: How to check cache version vs container and ref genome? Prefer user-supplied cache 116,
    //       as that is the version of the container.

    //
    // Run Ensembl VEP or rust-based fastpvep optionally
    //

    // Method is controlled by modules.config via ext.when, only one method is used
    ENSEMBLVEP_VEP(
        vcf_for_vep,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        fasta,
    )

    // FASTVEP is not implemented yet
    // FASTVEP_ANNOTATE()

}