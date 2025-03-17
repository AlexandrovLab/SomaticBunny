nextflow.enable.dsl=2

process RECALIBRATE_BQSR {
    scratch true
    label 'RECALIBRATE'
    conda "${params.java_env}"
    publishDir("${params.recal_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3
        

    input: 
    tuple val(patient), val(meta), val(status), path(bam), path(table), path(bai)
    each chunk


    output:
    tuple val(patient), val(meta.status), val(meta), path("*bam"), emit: MergeBam_input
    path("*bai"), emit: bai

    script:
    """
    /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    -L ${params.database_dir}/interval_list_20/${chunk}-scattered.interval_list \
    --bqsr-recal-file ${table} \
    -O ${meta.patient}_${meta.status}_recalibrated_${chunk}.bam
    """

}
