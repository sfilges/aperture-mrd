process LOFREQ_SOMATIC {
    tag "$meta.id"
    label 'process_high'

    container 'biocontainers/lofreq:2.1.5--py38h588ecb2_4'

    publishDir "${params.outdir}/${workflow.runName}/variant_calling/lofreq/{meta.id}", mode: params.publish_dir_mode, pattern: "*{vcf.gz,vcf.gz.tbi,log}"

    input:
    tuple val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Lofreq expected bam files, not cram
    def tumor_cram =  tumor.Extension == "cram" ? true : false
    def normal_cram =  normal.Extension == "cram" ? true : false

    def tumor_out = tumor_cram ? tumor.BaseName + ".bam" : "${tumor}"
    def normal_out = normal_cram ? normal.BaseName + ".bam" : "${normal}"

    // Builds the command to convert CRAM to BAM if necessary
    def samtools_cram_convert = ''
    samtools_cram_convert += normal_cram ? "    samtools view -T $fasta $normal -@ $task.cpus -o $normal_out\n" : ''
    samtools_cram_convert += normal_cram ? "    samtools index $normal_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools view -T ${fasta} $tumor -@ $task.cpus -o $tumor_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools index ${tumor_out}\n" : ''

    // Builds the command to remove CRAM files after variant calling
    def samtools_cram_remove = ''
    samtools_cram_remove += tumor_cram ? "     rm $tumor_out\n" : ''
    samtools_cram_remove += tumor_cram ? "     rm ${tumor_out}.bai\n " : ''
    samtools_cram_remove += normal_cram ? "     rm $normal_out\n" : ''
    samtools_cram_remove += normal_cram ? "     rm ${normal_out}.bai\n " : ''
    """
    $samtools_cram_convert

    lofreq \\
        somatic \\
        --threads $task.cpus \\
        $args \\
        -f $fasta \\
        -t $tumor_out \\
        -n $normal_out \\
        -d $dbsnp \\
        -o ${prefix}

    $samtools_cram_remove

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}