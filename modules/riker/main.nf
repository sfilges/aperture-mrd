process RIKER_MULTI {
    // https://github.com/fulcrumgenomics/riker
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/samtools:1.19.2--h50ea8bc_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai
    path  "versions.yml"           , emit: versions

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