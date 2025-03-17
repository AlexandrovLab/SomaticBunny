nextflow.enable.dsl=2

process RECALIBRATE_MergeReport {
    scratch true
    label 'process_low'
    conda "${params.java_env}"
    publishDir("${params.recal_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3
        

    input:
    val(map)

    output:
    tuple val(map.patient), val(map.meta), val(map.status), val(map.bam), path("*_merged_recal.table"), val(map.bai), emit: BQSR_input

    script:
    def cmd = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk GatherBQSRReports --tmp-dir . -O ${map.patient}_${map.status}_merged_recal.table"

    for( int i=0; i<20; i++ ) {
        cmd += " -I ${map.table[i]} "
    }

    cmd
}
