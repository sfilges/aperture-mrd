process GATK4_ADDORREPLACEREADGROUPS {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/gatk4:'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.rg.bam"), emit: rg_bam

    script:
    """
    gatk --java-options "-Xms4g -Xmx8g" \\
        AddOrReplaceReadGroups \\
        --INPUT $bam \\
        --OUTPUT ${meta.id}.rg.bam \\
        --RGLB ${meta.id} \\
        --RGPL illumina \\
        --RGPU ${meta.id} \\
        --RGSM ${meta.id} \\
        --SORT_ORDER coordinate \\
        --CREATE_INDEX true
    """

}