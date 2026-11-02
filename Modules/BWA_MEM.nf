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
    def batch_bases  = 160000000L                      // 16 CPUs x 10 Mb
    def rg_sm        = "${meta.patient}_${meta.sample}_${meta.status}"  
    """
    # PU = <flowcell>.<lane> of the first read (Illumina read names); falls back to the sample name
    pu=\$(gzip -cdf ${meta.fastq_1} | head -n 1 | awk -F'[:[:space:]]' '/^@/ && NF>=7 {print \$3"."\$4}')
    [ -n "\$pu" ] || pu=${rg_sm}

    set -o pipefail

    bwa-mem2 mem -Y -T 30 -t ${bwa_threads} -K ${batch_bases} \\
        -R "@RG\\tID:${rg_sm}\\tSM:${rg_sm}\\tLB:${rg_sm}\\tPL:ILLUMINA\\tPU:\$pu" \\
        ${params.ref} ${meta.fastq_1} ${meta.fastq_2} \\
    | samtools sort -@ ${sort_threads - 1} -T ${rg_sm} -o ${rg_sm}_raw.bam -

    """

}