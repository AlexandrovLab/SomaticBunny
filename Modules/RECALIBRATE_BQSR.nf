nextflow.enable.dsl=2

process RECALIBRATE_BQSR {
    conda "${params.java_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy = 'retry'
    maxRetries 3
        
    input: 
    tuple val(patient), val(meta), path(bam), path(table), path(bai)
    each chunk

    output:
    tuple val(patient), val(meta), path("*bam"), emit: MergeBam_input


    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    -L ${params.database_dir}/interval_list_20/${chunk}-scattered.interval_list \
    --bqsr-recal-file ${table} \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recalibrated_${chunk}.bam
    """

}
