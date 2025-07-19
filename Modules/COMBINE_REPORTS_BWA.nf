nextflow.enable.dsl=2

process COMBINE_REPORTS_BWA {
    label 'process_low'
    publishDir("${params.report_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3
    
    input:
    path individual_reports_bwa

    output:
    path "*.txt", emit: final_report_bwa

    script:
    """
    cat ${individual_reports_bwa.join(' ')} > quickcheck_bwa_report.txt
    """
}