nextflow.enable.dsl=2

process MERGE_HEADER {
    conda "${params.POST_env}"
    scratch true
    label 'process_low'
    publishDir("${params.post_dir}/standard_vcf/${patient_sample}/*_standard_PASSed.vcf")

    input:
    tuple val(map), val(scenario)

    output:
    path("*_standard_PASSed.vcf"), emit: final_standard_vcf
    script:
    """
    set -e 

    bash ${projectDir}/rename_info.sh \
        ${map[1]} \
        ${map[3]} \
        ${muse_vcf[4]}

    bash ${projectDir}/merge_header.sh \
        output_mutect.vcf \
        output_strelka.vcf \
        output_muse.vcf \
        ${sage_vcf} > merged_header.vcf

    cat merged_header.vcf ${final_vcf} > ${sample}_${scenario}_standard_PASSed.vcf
    """

}
