nextflow.enable.dsl=2

process MUTECT2_CALLING {
    scratch true
    conda "${params.gatk_env}"
    label 'MUTECT2_CALLING'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)
    each chunk
    path interval_dir

    output:
    tuple val(map), path("*f1r2.tar.gz"), emit: LearnReadOrientationModel_input
    tuple val(map), path("*.stats"), emit: MergeMutectStats_input
    tuple val(map), path("*vcf"), emit: MergeVcfs_input
    tuple val(map), path("*.idx")
    tuple val(map), path("*vcf"), path("*.idx"), path("*.stats"), path("*f1r2.tar.gz"), emit: mutect_call

    script:
    """
    ${params.gatk} Mutect2 -independent-mates -R ${params.ref} -pon ${params.mutect2_pon} --germline-resource ${params.mutect2_germline} --native-pair-hmm-threads ${task.cpus} --af-of-alleles-not-in-resource 0.00003125 --f1r2-tar-gz ${map.patient}_${map.tumor_meta.sample}_${chunk}-f1r2.tar.gz --normal-sample ${map.patient}_${map.normal_meta.sample}_normal --input ${map.normal} --tumor-sample ${map.patient}_${map.tumor_meta.sample}_tumor --input ${map.tumor} -L ${interval_dir}/${chunk}-scattered.interval_list -O ${map.patient}_${map.tumor_meta.sample}_${chunk}-unfiltered.vcf
    """
}
