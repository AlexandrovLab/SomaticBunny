nextflow.enable.dsl=2

process GETpileUP {
    conda "${params.java_env}"
    scratch true
    conda "${params.java_env}"
    label 'process_medium'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)
    each chunk

    output:
    tuple val(patient), val(meta), path("*.table"), val(chunk), emit: GETpileUP_Merge_input

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk GetPileupSummaries --java-options \"-Xmx\$(free -h | grep Mem | awk '{split(\$7,a,\"G\"); if(a[1]>5) print a[1]-5\"G\"; else print \"4G\"}')\" -I ${bam} -V ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz -L ${params.database_dir}/interval_list_20/${chunk}-scattered.interval_list -O ${meta.patient}_${meta.sample}_getpileupsummaries_${meta.status}_${chunk}.table
    """
    
}
