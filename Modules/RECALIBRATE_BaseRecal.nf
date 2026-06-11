nextflow.enable.dsl=2

process RECALIBRATE_BaseRecal {
    conda "${params.java_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy = 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)
    each chunk
    path interval_dir

    output:
    tuple val(patient), val(meta), path(bam), path("*table"), path(bai), emit: MergeReport_input

    script:
    """
    ${params.database_path}/gatk-4.6.0.0/gatk BaseRecalibrator \
    -I ${bam} \
    -R ${params.ref} \
    --known-sites ${params.recal_knownsite1} \
    --known-sites ${params.recal_knownsite2} \
    -L ${interval_dir}/${chunk}-scattered.interval_list \
    -O ${meta.patient}_${meta.sample}_${meta.status}_${chunk}_recal.table

    """

}
