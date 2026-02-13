process MuSE {
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
        ${params.MuSE2} call -f ${params.ref} -n 8 -O ${map.patient}_${map.tumor_meta.sample} ${map.tumor} ${map.normal}

        ${params.MuSE2} sump -I ${map.patient}_${map.tumor_meta.sample}.MuSE.txt -n 8  -E -O ${map.patient}_${map.tumor_meta.sample}.vcf ${dbsnp_param}
        """
    else
        """
        ${params.MuSE2} call -f ${params.ref} -n 8  -O ${map.patient}_${map.tumor_meta.sample} ${map.tumor} ${map.normal}

        ${params.MuSE2} sump -I ${map.patient}_${map.tumor_meta.sample}.MuSE.txt -n 8 -G -O ${map.patient}_${map.tumor_meta.sample}.vcf ${dbsnp_param}
        """
}

