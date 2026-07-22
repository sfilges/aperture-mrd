process MULTIQC {
    label 'process_single'

    container 'biocontainers/multiqc:1.35--pyhdfd78af_1'

    input:
    path(multiqc_files, stageAs: "?/*")

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data
    path "*_plots",              emit: plots, optional: true
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed -e "s/multiqc, version //g"'), emit: versions_multiqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    multiqc \\
        --force \\
        $args \\
        .
    """
}
