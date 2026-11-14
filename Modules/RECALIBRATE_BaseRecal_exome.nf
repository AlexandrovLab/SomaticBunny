nextflow.enable.dsl=2

process RECALIBRATE_BaseRecal_exome {
    conda "${params.gatk_env}"
    scratch true
    label 'RECALIBRATE'
    errorStrategy = 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)

    output:
    tuple val(patient), val(meta), path(bam), path("*table"), path(bai), emit: BQSR_input

    script:
    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    """
    ${params.gatk} BaseRecalibrator --java-options "-Xmx${xmx}M ${gc}" \
    -I ${bam} \
    -R ${params.ref} \
    --known-sites ${params.recal_knownsite1} \
    --known-sites ${params.recal_knownsite2} \
    -L ${params.recal_interval_wes} \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recal.table
    """
}
