nextflow.enable.dsl=2

process RECALIBRATE_BQSR {
    conda "${params.gatk_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy 'retry'
    maxRetries 3
        
    input: 
    tuple val(patient), val(meta), path(bam), path(table), path(bai)
    each chunk
    path interval_dir

    output:
    tuple val(patient), val(meta), path("${meta.patient}_${meta.sample}_${meta.status}_recalibrated_${chunk}.bam"), emit: MergeBam_input


    script:
    def unmapped = chunk == 20 ? '-L unmapped' : ''
    """
    ${params.gatk} ApplyBQSR \
    -R ${params.ref} \
    -I ${bam} \
    -L ${interval_dir}/bqsr/${chunk}-contigs.interval_list ${unmapped} \
    --bqsr-recal-file ${table} \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recalibrated_${chunk}.bam

    """

}
