process FULCRUM_FGUMI{
    // https://github.com/fulcrumgenomics/fgumi
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/fgumi:0.4.0--h54198d6_0'

}