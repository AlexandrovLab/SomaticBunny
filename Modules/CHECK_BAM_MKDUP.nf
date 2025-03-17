nextflow.enable.dsl=2

process CHECK_BAM_MKDUP {
    label 'process_low'
    conda "${params.samtools_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path("*.txt"), emit: individual_reports_mkdup

    script:
    """
    echo "Sample: ${meta.patient}_${meta.status}" >> ${meta.patient}_${meta.status}_quickcheck_mkdup.txt
    echo "Checking: ${bam}" >> ${meta.patient}_${meta.status}_quickcheck_mkdup.txt
    result=\$(samtools quickcheck ${bam} 2>&1)
    if [ -z "\$result" ]; then
        echo "PASS: No issues detected" >> ${meta.patient}_${meta.status}_quickcheck_mkdup.txt
    else
        echo "FAIL: Issues detected" >> ${meta.patient}_${meta.status}_quickcheck_mkdup.txt
        echo "\$result" >> ${meta.patient}_${meta.status}_quickcheck_mkdup.txt
    fi
    echo "------------------------" >> ${meta.patient}_${meta.status}_quickcheck_mkdup.txt
    """
}