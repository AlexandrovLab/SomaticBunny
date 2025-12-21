nextflow.enable.dsl=2

process POST {
    scratch true
    label 'process_low'  
    publishDir("${params.post_dir}/${patient_sample}", mode: 'copy', pattern: '*_PASSed.vcf')
    
    input:
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai)

    output:
    path("snvs_filtered/2outof4/${patient_sample}_PASSed.vcf"), emit: snv_passed
    path("indels_filtered/2outof3/${patient_sample}_PASSed.vcf"), emit: indel_passed

    script:
    """
    # Create symlink so pysam can find the index file
    ln -sf ${tumor_bai} ${tumor_bam}.bai

    bash ${params.database_path}/EVC_nextflow/Databases/filtering.sh \
        ${patient_sample} \
        ${mutect2_vcf} \
        ${muse_vcf} \
        ${strelka_snv} \
        ${strelka_indel} \
        ${sage_vcf} \
        ${tumor_bam} \
        ${params.database_path} \
        ${params.ref}
    """
}
