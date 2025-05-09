nextflow.enable.dsl=2

process RECALIBRATE_SortBam {
    scratch true
    label 'RECALIBRATE_SortBam'
    conda "${params.samtools_env}"
    if (params.publish.tokenize(',').contains('recal_bam')) {
        publishDir("${params.recal_dir}", mode: 'copy')
    }
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam)

    output:
    tuple val(patient), val(meta), path("*bam"), path("*bai"), emit: pair_recal
    tuple val(meta), path("*bam"), path("*bai"), emit: recal_bam

    script:
    def sorted_bam = "${meta.patient}_${meta.sample}_${meta.status}_recal.bam"
    """
    samtools sort -@ 4 ${bam} -o ${sorted_bam}
    samtools index ${sorted_bam}
    """
}



