nextflow.enable.dsl=2

process STRELKA {
    conda "${params.strelka_env}"
    scratch true
    label 'STRELKA'
    publishDir("${params.strelka_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    path("*.somatic_indels.vcf.gz"), emit: STRELKA_indels
    path("*.somatic_indels.vcf.gz.tbi"), emit: STRELKA_indels_tbi
    path("*.somatic_snvs.vcf.gz"), emit: STRELKA_snvs
    path("*.somatic_snvs.vcf.gz.tbi"), emit: STRELKA_snvs_tbi
    val(map), emit: STRELKA_out


    script:
    if (map.type == "exome")
        """
        configureStrelkaSomaticWorkflow.py --referenceFasta ${params.ref} --normalBam ${map.normal} --tumorBam ${map.tumor} --runDir strelka --exome

        python2 strelka/runWorkflow.py -m local -j 10

        mv strelka/results/variants/somatic.indels.vcf.gz     ${map.patient}.somatic_indels.vcf.gz
        mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${map.patient}.somatic_indels.vcf.gz.tbi
        mv strelka/results/variants/somatic.snvs.vcf.gz       ${map.patient}.somatic_snvs.vcf.gz
        mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${map.patient}.somatic_snvs.vcf.gz.tbi
        """
    else
        """
        configureStrelkaSomaticWorkflow.py --referenceFasta ${params.ref} --normalBam ${map.normal} --tumorBam ${map.tumor} --runDir strelka --callRegions ${params.bed}

        python2 strelka/runWorkflow.py -m local -j 10

        mv strelka/results/variants/somatic.indels.vcf.gz     ${map.patient}.somatic_indels.vcf.gz
        mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${map.patient}.somatic_indels.vcf.gz.tbi
        mv strelka/results/variants/somatic.snvs.vcf.gz       ${map.patient}.somatic_snvs.vcf.gz
        mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${map.patient}.somatic_snvs.vcf.gz.tbi
        """
}