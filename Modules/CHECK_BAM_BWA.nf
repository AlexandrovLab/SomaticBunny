nextflow.enable.dsl=2

/*
 * Runs samtools quickcheck on one BAM and reports PASS or FAIL.
 *
 * The task itself always succeeds when quickcheck runs, so every sample gets
 * checked and gets a report. main.nf passes only PASS BAMs on to MarkDuplicates.
 *
 * Outputs:
 *   checked                 tuple(meta, bam, status)   status = "PASS" or "FAIL"
 *   individual_reports_bwa  <patient>_<sample>_<status>_quickcheck_bwa.txt
 */
process CHECK_BAM_BWA {
    label 'process_low'
    conda "${params.samtools_env}"
    // Retry only for infrastructure problems (killed job, out of memory, etc.).
    // A corrupt BAM is reported as FAIL, not as a task error, so it is never retried.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), env(QC_STATUS), emit: checked
    path("*_quickcheck_bwa.txt"),               emit: individual_reports_bwa

    script:
    def report = "${meta.patient}_${meta.sample}_${meta.status}_quickcheck_bwa.txt"
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
