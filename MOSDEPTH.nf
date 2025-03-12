nextflow.enable.dsl=2

process MOSDEPTH {
    conda "${params.mosdepth_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.mosdepth_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3
        
    input:
    val(map)

    output:
    tuple val(map.patient), path("*txt"), emit: coverage
    tuple val(map.patient), path("*txt"), emit: MOSDEPTH_out

    script:
    if (map.type == "exome")
        """
        mosdepth -t 8 --by ${params.database_dir}/GRCh38_exome.bed ${map.patient}_tumor ${map.tumor}
        mosdepth --by ${params.database_dir}/GRCh38_exome.bed ${map.patient}_normal ${map.normal}
        """

    else
        """
        mosdepth -t 8 -n --fast-mode --by 500 ${map.patient}_tumor ${map.tumor}
        mosdepth -n --fast-mode --by 500 ${map.patient}_normal ${map.normal}
        """
}