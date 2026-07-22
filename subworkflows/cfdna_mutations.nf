/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CFDNA_MUTATIONS Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Process WES/WGS data from cfDNA in tumor-informed or tumor-agnostic modes.

    tumor-informed:
        - cfDNA cram
        - VCF of the tumor-specfic SNV compendium

    tumor-agnostic:
        - cfDNA cram

    In WGS mode we make use of a read-level error model and genome-wide integration for SNVs. 
    WES mode performs traditional 'locus-centric' calling.

    Perform CNV calling:
        - ichorCNA on WGS to estimate tumor fraction (cfDNA tumor fractio  >1-3%)
        - cnvkit on WES

    MSI:
        - msisensor2 (cfDNA tumor fractio  >0.1%)

    Fragmentomic features are processed using a separate workflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CFDNA_MUTATIONS {
    // stub
}