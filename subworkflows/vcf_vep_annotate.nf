/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT ANNOTATION Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Annotates VCF files using Enembl VEP.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ENSEMBLVEP_VEP } from '../modules/annotation/ensemblvep/vep/main'

workflow VCF_VEP_ANNOTATE {
    take:
        vcf           // tuple val(meta), path(vcf)
        genome        // params.vep_genome
        species       // params.vep_species
        cache_version // params.vep_cache_version
        cache         // params.vep_cache
        fasta         // tuple val(meta2), path(fasta)

    main:
    // TODO: Need to download cache first if not available (in separate workflow?)
    // TODO: How to check cache version vs container and ref genome? Prefer user-supplied cache 116,
    //       as that is the version of the container.

    //
    // Run Ensembl VEP or rust-based fastpvep optionally
    //
    ENSEMBLVEP_VEP(
        vcf,
        genome,
        species,
        cache_version,
        cache,
        fasta,
    )

}