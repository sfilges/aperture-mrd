process MULTIQC {
    label 'process_single'

    container 'biocontainers/multiqc:1.25.2--pyhdfd78af_0'

    publishDir "${params.outdir}/${workflow.runName}/reports/multiqc", mode: params.publish_dir_mode

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
