nextflow.enable.dsl=2

process GETpileUP {
    scratch true
    conda "${params.gatk_env}"
    label 'process_medium'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)
    each chunk
    path interval_dir

    output:
    tuple val(patient), val(meta), path("*.table"), val(chunk), emit: GETpileUP_Merge_input

    script:
    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    """
    ${params.gatk} GetPileupSummaries --java-options "-Xmx${xmx}M" -I ${bam} -V ${params.mutect2_germline} -L ${interval_dir}/${chunk}-scattered.interval_list -O ${meta.patient}_${meta.sample}_getpileupsummaries_${meta.status}_${chunk}.table
    """
}


