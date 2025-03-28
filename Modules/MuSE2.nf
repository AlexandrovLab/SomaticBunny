process MuSE2 {
    conda "${params.muse2_env}"
    publishDir("${params.muse2_dir}", mode: 'copy')
    conda "${params.muse_env}"
    scratch true
    label 'MuSE2'
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    path("*.vcf"), emit: MuSE2_vcf
    val(map), emit:MuSE2_out

    script:
    if (params.type == "exome")
        """
        ${params.MuSE2} call -f ${params.ref} -n 8 -O ${map.patient} ${map.tumor} ${map.normal}

        ${params.MuSE2} sump -I ${map.patient}.MuSE.txt -n 8 -E -O ${map.patient}_${map.tumor_meta.sample}.vcf -D ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz
        """
    else
        """
        ${params.MuSE2} call -f ${params.ref} -n 16 -O ${map.patient} ${map.tumor} ${map.normal}

        ${params.MuSE2} sump -I ${map.patient}.MuSE.txt -n 16 -G -O ${map.patient}.vcf -D ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz
        """
}
