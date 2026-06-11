nextflow.enable.dsl=2

process CalculateContamination {
    scratch true
    label 'process_low'
    conda "${params.java_env}"
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map), path("*contamination.table"), emit: MUTECT2_contamination_table
    tuple val(map), path("*segments.table"), emit: MUTECT2_segments_table

    script:
    """
    ${params.database_path}/gatk-4.6.0.0/gatk CalculateContamination -I ${map.tumor} -matched ${map.normal} -O ${map.patient}_${map.tumor_meta.sample}_contamination.table --tumor-segmentation ${map.patient}_${map.tumor_meta.sample}_segments.table
    """
}
