process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0'

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val   genome
    val   species
    val   cache_version
    path  cache
    tuple val(meta2), path(fasta)
    path  extra_files

    output:
    tuple val(meta), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(meta), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(meta), path("*.json.gz") , optional:true, emit: json
    path "*.summary.html"              , optional:true, emit: report
    tuple val("${task.process}"), val('ensemblvep'), eval("echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//'"), emit: versions_ensemblvep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.tab.gz
    touch ${prefix}.json.gz
    touch ${prefix}.summary.html
    """
}