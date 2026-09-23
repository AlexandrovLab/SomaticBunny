nextflow.enable.dsl=2

process SUMMARY{
    conda "${params.summary_env}"
    publishDir("$projectDir", mode: 'copy')
    label 'process_low'
    errorStrategy = 'retry'
    maxRetries 3
    input:
    path saved_csv
    path sample_csv

    output:
    stdout

    script:
    """
    python ${projectDir}/summary.py "${saved_csv}" "${sample_csv}"
    """

}