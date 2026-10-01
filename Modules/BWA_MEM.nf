nextflow.enable.dsl=2

process BWA_MEM {
    conda "${params.bwamem2_env}"
    scratch true
    label 'BWA_MEM'
    if (params.publish.tokenize(',').contains('raw_bam')) {
        publishDir("${params.bam_dir}", mode: 'copy')
    }
    maxRetries 3

    input:
    val(meta)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def sort_threads = Math.max(1, task.cpus.intdiv(4))         // a quarter of task.cpus, at least 1
    def bwa_threads  = Math.max(1, task.cpus - sort_threads)     // the rest, at least 1
    def batch_bases  = 10000000L * task.cpus                     // same batches as bwa-mem2 -t task.cpus
    """
    set -o pipefail

    bwa-mem2 mem -T 0 -t ${bwa_threads} -K ${batch_bases} -R "@RG\\tID:${meta.patient}\\tSM:${meta.patient}_${meta.sample}_${meta.status}\\tPL:ILLUMINA" ${params.ref} ${meta.fastq_1} ${meta.fastq_2} \
    | samtools sort -@ ${sort_threads - 1} -T ${meta.patient}_${meta.sample}_${meta.status} -o ${meta.patient}_${meta.sample}_${meta.status}_raw.bam -

    """
}