nextflow.enable.dsl=2

process RECALIBRATE_BQSR_exome {
    conda "${params.gatk_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy = 'retry'
    maxRetries 3
        
    input: 
    tuple val(patient), val(meta), path(bam), path(table), path(bai)

    output:
    tuple val(patient), val(meta),
          path("*_recal.bam"),
          path("*.bai"),
          emit: pair_recal

    script:
    """
    ${params.gatk} ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    --bqsr-recal-file ${table} \
    --create-output-bam-index true \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recal.bam
    """
}
