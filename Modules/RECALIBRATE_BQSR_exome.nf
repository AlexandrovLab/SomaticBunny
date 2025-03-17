nextflow.enable.dsl=2

process RECALIBRATE_BQSR_exome {
    scratch true
    label 'RECALIBRATE'
    conda "${params.java_env}"
    publishDir("${params.recal_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3
        

    input: 
    tuple val(patient), val(meta), val(status), path(bam), path(table), path(bai)


    output:
    tuple val(patient), val(meta.status), val(meta), path("*bam"), emit: SortBam_input
    path("*bai"), emit: bai

    script:
    """
    /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    -L ${params.database_dir}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list \
    --bqsr-recal-file ${table} \
    -O ${meta.patient}_${meta.status}_recalibrated.bam
    """

}
