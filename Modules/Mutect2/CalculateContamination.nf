nextflow.enable.dsl=2

process CalculateContamination {
    conda "${params.java_env}"
    scratch true
    label 'process_low'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map.patient), path("*contamination.table"), emit: MUTECT2_contamination_table
    tuple val(map.patient), path("*segments.table"), emit: MUTECT2_segments_table

    script:
    """
    /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk CalculateContamination -I ${map.tumor} -matched ${map.normal} -O ${map.patient}_contamination.table --tumor-segmentation ${map.patient}_segments.table
    """
}
