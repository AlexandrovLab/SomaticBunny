nextflow.enable.dsl=2

process SAGE {
    conda "${params.SAGE_java_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.SAGE_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    path("*vcf.gz"), emit: SAGE_vcf_gz
    path("*vcf.gz.tbi"), emit: SAGE_vcf_gz_tbi
    tuple val(map),path("*vcf.gz"), emit: SAGE_out

    script:
    """
    java -Xms4G -Xmx32G -cp ${params.SAGE_jar} com.hartwig.hmftools.sage.SageApplication \
        -threads ${task.cpus} \
        -reference ${map.patient}_${map.normal_meta.sample}_normal \
        -reference_bam  ${map.normal}\
        -tumor ${map.patient}_${map.tumor_meta.sample}_tumor \
        -tumor_bam ${map.tumor} \
        -ref_genome_version ${params.SAGE_ref_genome_version} \
        -ref_genome ${params.ref} \
        -hotspots ${params.SAGE_hotspots} \
        -panel_bed ${params.SAGE_panel_bed} \
        -high_confidence_bed ${params.SAGE_high_confidence_bed} \
        -ensembl_data_dir ${params.SAGE_ref_dir}/common/ensembl_data \
        -out ${map.patient}_${map.tumor_meta.sample}.sage.vcf.gz
    """
}

