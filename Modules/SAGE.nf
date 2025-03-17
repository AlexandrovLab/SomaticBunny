nextflow.enable.dsl=2

process SAGE {
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
    val(map), emit: SAGE_out

    script:
    """
    java -Xms4G -Xmx32G -cp ${params.SAGE_java} com.hartwig.hmftools.sage.SageApplication \
        -threads 8 \
        -reference ${map.patient}_normal \
        -reference_bam   ${map.normal}\
        -tumor ${map.patient}_tumor \
        -tumor_bam ${map.tumor} \
        -ref_genome_version 38 \
        -ref_genome ${params.ref} \
        -hotspots ${params.SAGE_ref_dir}/variants/KnownHotspots.somatic.38.vcf.gz \
        -panel_bed ${params.SAGE_ref_dir}/variants/ActionableCodingPanel.38.bed.gz \
        -high_confidence_bed ${params.SAGE_ref_dir}/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz \
        -ensembl_data_dir ${params.SAGE_ref_dir}/common/ensembl_data \
        -out ${map.patient}.sage.vcf.gz
    """

}