nextflow.enable.dsl=2

process MKDUP {
    conda "${params.mkdup_env}"
    scratch true
    label 'MKDUP'
    publishDir("${params.mkdup_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3
        
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta.patient), val(meta), path("*bam"), path("*bai"), emit: pair_mutect
    tuple val(meta), path("*bam"), path("*bai"), emit: mkdup_bam

    script:
    """
    picard MarkDuplicates ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES=4000 MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true -Djava.io.tmpdir=${params.mkdup_temp_dir} -XX:ParallelGCThreads=8 -Xmx16g VALIDATION_STRINGENCY=STRICT I=${bam} O=${meta.patient}_${meta.status}_mkdp.bam M=${meta.patient}_${meta.status}_markDuplicates_Matrix.txt
    """
}