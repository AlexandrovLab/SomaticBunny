nextflow.enable.dsl=2

process CHECK_BAM_BWA {
    label 'process_low'
    conda "${params.samtools_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), path(bam)

    output:
    path("*.txt"), emit: individual_reports_bwa

    script:
    """
    echo "Sample: ${meta.patient}_${meta.status}" >> ${meta.patient}_${meta.status}_quickcheck_bwa.txt
    echo "Checking: ${bam}" >> ${meta.patient}_${meta.status}_quickcheck_bwa.txt
    result=\$(samtools quickcheck ${bam} 2>&1)
    if [ -z "\$result" ]; then
        echo "PASS: No issues detected" >> ${meta.patient}_${meta.status}_quickcheck_bwa.txt
    else
        echo "FAIL: Issues detected" >> ${meta.patient}_${meta.status}_quickcheck_bwa.txt
        echo "\$result" >> ${meta.patient}_${meta.status}_quickcheck_bwa.txt
    fi
    echo "------------------------" >> ${meta.patient}_${meta.status}_quickcheck_bwa.txt
    """
}