nextflow.enable.dsl=2

process COMBINE_REPORTS_RECAL {
    label 'process_low'
    publishDir("${params.report_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path individual_reports_recal

    output:
    path "*.txt", emit: final_report_recal

    script:
    """
    cat ${individual_reports_recal.join(' ')} > quickcheck_recal_report.txt
    """
}