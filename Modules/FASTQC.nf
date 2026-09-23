nextflow.enable.dsl=2

process FASTQC {
    conda "${params.fastqc_env}"
    scratch true
    label 'FASTQC'
    publishDir("${params.FASTQC_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(meta)

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    tuple val(meta), path("*.html"), emit: fastqc_out

    script:
    """
    fastqc -t ${task.cpus} -o ./ ${meta.fastq_1} ${meta.fastq_2}
    """
}