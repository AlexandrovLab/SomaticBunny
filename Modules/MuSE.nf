process MuSE {
    conda "${params.muse_env}"
    publishDir("${params.muse_dir}", mode: 'copy')
    scratch true
    label 'MuSE'
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    path("*.vcf"), emit: MuSE_vcf
    tuple val(map), path("*.vcf"), emit:MuSE_out

    script:
    def dbsnp_param = params.muse_vcf ? "-D ${params.muse_vcf}" : ""
    
    if (params.type == "exome")
        """
        MuSE call -f ${params.ref} -O ${map.patient}_${map.tumor_meta.sample} ${map.tumor} ${map.normal}

        MuSE sump -I ${map.patient}_${map.tumor_meta.sample}.MuSE.txt -E -O ${map.patient}_${map.tumor_meta.sample}.vcf ${dbsnp_param}
        """
    else
        """
        MuSE call -f ${params.ref} -O ${map.patient}_${map.tumor_meta.sample} ${map.tumor} ${map.normal}

        MuSE sump -I ${map.patient}_${map.tumor_meta.sample}.MuSE.txt -G -O ${map.patient}_${map.tumor_meta.sample}.vcf ${dbsnp_param}
        """
}

