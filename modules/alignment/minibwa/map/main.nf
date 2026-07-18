process MINIBWA_MAP {
    // https://github.com/lh3/minibwa
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/minibwa:0.4--hab16a5f_0'
}