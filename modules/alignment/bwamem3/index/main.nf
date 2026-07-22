process BWAMEM3_INDEX {
    // https://bwa-mem3.readthedocs.io/en/v0.6.0/
    tag "${meta.id}"

    // No process_high label: its withLabel memory (72.GB) would win over this in-script
    // directive. cpus/time are set via withName in conf/modules.config instead.
    // 64-bit SA over fwd+rev of hg38 peaks ~70 GiB; scale memory with genome size.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    container 'community.wave.seqera.io/library/bwa-mem3_htslib_samtools:391ed2ac52c4a15a'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem3"), emit: index
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // ~85% of the allocation: bwa-mem3's default cap (32 GiB) aborts the build otherwise
    def max_mem = Math.max(1, Math.round(task.memory.toGiga() * 0.85))
    """
    mkdir bwamem3

    bwa-mem3 \\
        index \\
        -p bwamem3/${prefix} \\
        -t ${task.cpus} \\
        --max-memory ${max_mem}g \\
        ${args} \\
        ${fasta}
    """
}