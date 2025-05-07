nextflow.enable.dsl=2

process FilterMutectCalls {
    scratch true
    label 'process_low'
    conda "${params.java_env}"
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy = 'retry'
    maxRetries 3

    input:
    tuple val(map),
          path(unfiltered_vcf),
          path(contamination_table),
          path(segments_table),
          path(read_orientation_model_tar),
          path(merged_stats)

    output:
    tuple val(map), path("*vcf"), path("*idx"), path("*stats"), emit: MUTECT2_final_out
    tuple val (map),path("*vcf") emit: Mutect2_out

    script:
    """
    ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk FilterMutectCalls \
    -R ${params.ref} \
    -V ${unfiltered_vcf} \
    --contamination-table ${contamination_table} \
    --ob-priors ${read_orientation_model_tar} \
    -O ${map.patient}_${map.tumor_meta.sample}_mutect2_filtered.vcf \
    --stats ${merged_stats} \
    --filtering-stats ${map.patient}_${map.tumor_meta.sample}_mutect2_filtered.stats \
    --tumor-segmentation ${segments_table}
    """
}
