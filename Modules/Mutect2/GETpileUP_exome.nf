nextflow.enable.dsl=2

process GETpileUP_exome {
    scratch true
    conda "${params.gatk_env}"
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
    ${params.gatk} GetPileupSummaries --java-options \"-Xmx\$(free -h | grep Mem | awk '{split(\$7,a,\"G\"); if(a[1]>5) print a[1]-5\"G\"; else print \"4G\"}')\" -I ${bam} -V ${params.mutect2_germline} -L ${params.recal_interval_wes} -O ${meta.patient}_${meta.sample}_getpileupsummaries_${meta.status}.table
    """
    
}
