process FASTVEP {
    // https://github.com/Huang-lab/fastVEP
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/riker:0.4.1--hec9b1f2_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    riker multi \\
        -i sample.bam \\
        -r ref.fa \\
        -o out_prefix \\
        --tools alignment isize basic wgs gcbias alignment error \\
        ${args}
    """
}