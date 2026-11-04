nextflow.enable.dsl=2

/*
 * Runs samtools quickcheck on one BAM and reports PASS or FAIL.
 *
 * The task itself always succeeds when quickcheck runs, so every sample gets
 * checked and gets a report. main.nf passes only PASS BAMs on, then drops any patient left
 * without a passing normal and tumor before variant calling.
 *
 * Outputs:
 *   checked                   tuple(patient, meta, bam, bai, status)   status = "PASS" or "FAIL"
 *   individual_reports_recal  <patient>_<sample>_<status>_quickcheck_recal.txt
 */
process CHECK_BAM_RECAL {
    label 'process_low'
    conda "${params.samtools_env}"
    // Retry only for infrastructure problems (killed job, out of memory, etc.).
    // A corrupt BAM is reported as FAIL, not as a task error, so it is never retried.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)

    output:
    tuple val(patient), val(meta), path(bam), path(bai), env(QC_STATUS), emit: checked
    path("*_quickcheck_recal.txt"),                                      emit: individual_reports_recal

    script:
    def report = "${meta.patient}_${meta.sample}_${meta.status}_quickcheck_recal.txt"
    """
    {
        echo "Sample: ${meta.patient}_${meta.sample}_${meta.status}"
        echo "Checking: ${bam}"
    } > ${report}

    # quickcheck signals failure through its exit code; -v prints the failing file name.
    if samtools quickcheck -v ${bam} > quickcheck.log 2>&1; then
        QC_STATUS=PASS
        echo "PASS: No issues detected" >> ${report}
    else
        QC_STATUS=FAIL
        echo "FAIL: Issues detected" >> ${report}
        cat quickcheck.log >> ${report}
    fi
    echo "------------------------" >> ${report}
    """
}
