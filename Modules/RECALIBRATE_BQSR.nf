nextflow.enable.dsl=2

process RECALIBRATE_BQSR {
    conda "${params.java_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy 'terminate'
    maxRetries 3
        
    input: 
    tuple val(patient), val(meta), path(bam), path(table), path(bai)
    each chunk
    path interval_dir

    output:
    tuple val(patient), val(meta), path("*bam"), emit: MergeBam_input


    script:
    """
    ${params.database_path}/gatk-4.6.0.0/gatk ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    -L ${interval_dir}/${chunk}-scattered.interval_list \
    --bqsr-recal-file ${table} \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recalibrated_${chunk}.bam

    """

}
