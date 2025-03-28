nextflow.enable.dsl=2

process RECALIBRATE_BaseRecal_exome {
    conda "${params.java_env}"
    scratch true
    label 'RECALIBRATE'
    conda "${params.java_env}"
    publishDir("${params.recal_dir}", mode: 'copy')
    errorStrategy = 'retry'
    maxRetries 3
        
    input:
    tuple val(patient), val(meta), path(bam), path(bai)

    output:
    tuple val(patient), val(meta), path(bam), path("*table"), path(bai), emit: BQSR_input

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk BaseRecalibrator \
    -I ${bam} \
    -R ${params.ref} \
    --known-sites ${params.database_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ${params.database_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -L ${params.database_dir}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list \
    -O ${meta.patient}_${meta.sample}_${meta.status}_recal.table
    """
}
