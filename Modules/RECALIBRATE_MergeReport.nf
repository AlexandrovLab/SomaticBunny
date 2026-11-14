nextflow.enable.dsl=2

process RECALIBRATE_MergeReport {
    conda "${params.gatk_env}"
    scratch true
    label 'process_low'
    if (params.publish.tokenize(',').contains('recal_bam')) {
        publishDir("${params.recal_dir}", mode: 'copy')
    }
    errorStrategy = 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map.patient), val(map.meta), val(map.bam), path("*_merged_recal.table"), val(map.bai), emit: BQSR_input

    script:
    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    def cmd = "${params.gatk} GatherBQSRReports --java-options \"-Xmx${xmx}M ${gc}\" --tmp-dir . -O ${map.meta.patient}_${map.meta.sample}_${map.meta.status}_merged_recal.table"

    for( int i=0; i<20; i++ ) {
        cmd += " -I ${map.table[i]} "
    }

    cmd
}
