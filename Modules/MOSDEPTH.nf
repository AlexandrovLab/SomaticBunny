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
    tuple val(map.patient), path("*regions.bed.gz"), emit:regions_bed

    script:
    if (map.type == "exome")
        """
        mosdepth -t 8 --by ${params.database_dir}/GRCh38_exome.bed ${map.patient}_${map.tumor_meta.sample}_tumor ${map.tumor}
        mosdepth -t 8 --by ${params.database_dir}/GRCh38_exome.bed ${map.patient}_${map.tumor_meta.sample}_normal ${map.normal}
        """

    else
        """
        mosdepth -t 8 -n --fast-mode --by 500 ${map.patient}_${map.tumor_meta.sample}_tumor ${map.tumor}
        mosdepth -t 8  -n --fast-mode --by 500 ${map.patient}_${map.tumor_meta.sample}_normal ${map.normal}
        """
}
