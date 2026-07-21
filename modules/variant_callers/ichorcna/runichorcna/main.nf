process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'quay.io/dincalcilab/ichorcna:0.4.0-2ab0be2"'

    input:
    tuple val(meta), path(wig)
    path gc_wig
    path map_wig
    path panel_of_normals
    path centromere

    output:
    tuple val(meta), path("*.cna.seg")    , emit: cna_seg
    tuple val(meta), path("*.params.txt") , emit: ichorcna_params
    path "*genomeWide.pdf"                , emit: genome_plot
    tuple val("${task.process}"), val('ichorcna'), val('0.3.2'), emit: versions_ichorcna, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    """
    runIchorCNA.R \\
        $args \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir .

    cp */*genomeWide.pdf .
    """
}