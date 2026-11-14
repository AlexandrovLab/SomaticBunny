nextflow.enable.dsl=2

process CalculateContamination {
    scratch true
    label 'process_low'
    conda "${params.gatk_env}"
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map), path("*contamination.table"), emit: MUTECT2_contamination_table
    tuple val(map), path("*segments.table"), emit: MUTECT2_segments_table

    script:
    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    """
    ${params.gatk} CalculateContamination --java-options "-Xmx${xmx}M ${gc}" -I ${map.tumor} -matched ${map.normal} -O ${map.patient}_${map.tumor_meta.sample}_contamination.table --tumor-segmentation ${map.patient}_${map.tumor_meta.sample}_segments.table
    """
}
