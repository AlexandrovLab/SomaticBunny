nextflow.enable.dsl=2

process RECALIBRATE_BaseRecal {
    conda "${params.gatk_env}"
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
    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    """
    ${params.gatk} BaseRecalibrator --java-options "-Xmx${xmx}M ${gc}" \
    -I ${bam} \
    -R ${params.ref} \
    --known-sites ${params.recal_knownsite1} \
    --known-sites ${params.recal_knownsite2} \
    -L ${interval_dir}/${chunk}-scattered.interval_list \
    -O ${meta.patient}_${meta.sample}_${meta.status}_${chunk}_recal.table

    """

}
