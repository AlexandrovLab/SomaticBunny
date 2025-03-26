nextflow.enable.dsl=2

process MUTECT2_CALLING_exome {
    scratch true
    conda "${params.java_env}"
    label 'process_low'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map), path("*f1r2.tar.gz"), emit: LearnReadOrientationModel_input
    tuple val(map), path("*.stats"), emit: MUTECT2_stats
    tuple val(map), path("*vcf"), emit: MUTECT2_vcf
    tuple val(map), path("*.idx"), emit: MUTECT2_idx
    tuple val(map), path("*vcf"), path("*.idx"), path("*.stats"), path("*f1r2.tar.gz"), emit: mutect_call

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk Mutect2 -independent-mates -R ${params.ref} -pon ${params.database_dir}/MuTect2.PON.5210.vcf.gz --germline-resource ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz --native-pair-hmm-threads 2 --af-of-alleles-not-in-resource 0.0000025 --f1r2-tar-gz ${map.patient}_${map.tumor_meta.sample}-f1r2.tar.gz --normal-sample ${map.patient}_${map.tumor_meta.sample}_normal --input ${map.normal} --tumor-sample ${map.patient}_${map.tumor_meta.sample}_tumor --input ${map.tumor} -L ${params.database_dir}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list -O ${map.patient}_${map.tumor_meta.sample}-unfiltered.vcf
    """

}