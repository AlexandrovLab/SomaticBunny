nextflow.enable.dsl=2

process FASTQC {
    conda "${params.fastqc_env}"
    scratch true
    label 'FASTQC'
    publishDir("${params.FASTQC_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    tuple val(meta), path("*.html"), emit: fastqc_out

    script:
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    """
    export JAVA_TOOL_OPTIONS="${gc}"
    fastqc -t ${task.cpus} -o ./ ${fastq_1} ${fastq_2}
    """
}