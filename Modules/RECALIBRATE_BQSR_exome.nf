nextflow.enable.dsl=2

process RECALIBRATE_BQSR_exome {
    conda "${params.java_env}"
    scratch true
    label 'RECALIBRATE'
    if (params.publish.tokenize(',').contains('recal_bam')) {
        publishDir("${params.recal_dir}", mode: 'copy')
    }
    errorStrategy = 'retry'
    maxRetries 3
        
    input: 
    tuple val(patient), val(meta), path(bam), path(table), path(bai)

    output:
    tuple val(patient), val(meta), path("*bam"), emit: SortBam_input

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    -L ${params.database_dir}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list \
    --bqsr-recal-file ${table} \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recalibrated.bam
    """

}
