nextflow.enable.dsl=2

process GETpileUP_exome {
    conda "${params.java_env}"
    scratch true
    conda "${params.java_env}"
    label 'process_medium'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)

    output:
    tuple val(patient), val(meta), path("*.table"), emit: CalculateContamination_input

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk GetPileupSummaries --java-options \"-Xmx\$(free -h | grep Mem | awk '{split(\$7,a,\"G\"); if(a[1]>5) print a[1]-5\"G\"; else print \"4G\"}')\" -I ${bam} -V ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz -L ${params.database_dir}/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list -O ${meta.patient}_${meta.sample}_getpileupsummaries_${meta.status}.table
    """
    
}
