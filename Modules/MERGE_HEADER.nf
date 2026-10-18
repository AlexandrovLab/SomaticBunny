nextflow.enable.dsl=2

process MERGE_HEADER {
    conda "${params.POST_rename_env}"
    scratch true
    label 'process_low'
    publishDir("${params.post_dir}/standard_vcf/${map[0]}",mode:'copy',pattern: "*_standard_PASSed.vcf")

    input:
    tuple val(map), val(scenario)

    output:
    path("*_standard_PASSed.vcf"), emit: final_standard_vcf
    script:
    """
    bash ${projectDir}/merge_header.sh \
        ${map[1]} \
	${map[3]} \
        ${map[2]} \
        ${map[4]} > merged_header.vcf

    cat merged_header.vcf ${map[5]} > ${map[0]}_${scenario}_standard_PASSed.vcf
    """

}
