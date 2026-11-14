nextflow.enable.dsl=2

process COMBINE_REPORTS_MKDUP {
    label 'process_low'
    publishDir("${params.report_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3
    
    input:
    path individual_reports_mkdup

    output:
    path "*.txt", emit: final_report_mkdup

    script:
    """
    cat ${individual_reports_mkdup.join(' ')} > quickcheck_mkdup_report.txt
    """
}