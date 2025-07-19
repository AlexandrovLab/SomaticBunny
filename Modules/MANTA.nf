nextflow.enable.dsl=2

process MANTA {
    conda "${params.manta_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.manta_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map), path("*diploidSV.vcf.gz"), emit: MANTA_diploidSV
    tuple val(map), path("*somaticSV.vcf.gz"), emit: MANTA_somaticSV
    tuple val(map), path("*candidateSV.vcf.gz"), emit: MANTA_candidateSV
    tuple val(map), path("*candidateSmallIndels.vcf.gz"), emit: MANTA_candidateSmallIndels

    script:
    if (params.type == "exome")

        """
        configManta.py \
        --normalBam ${map.normal} \
        --tumorBam ${map.tumor} \
        --referenceFasta ${params.ref} \
        --runDir manta \
        --exome

        python2 manta/runWorkflow.py -j 8

        mv manta/results/variants/candidateSmallIndels.vcf.gz ${map.patient}_${map.sample}.candidateSmallIndels.vcf.gz
        mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi ${map.patient}_${map.sample}.candidateSmallIndels.vcf.gz.tbi

        mv manta/results/variants/candidateSV.vcf.gz ${map.patient}_${map.sample}.candidateSV.vcf.gz
        mv manta/results/variants/candidateSV.vcf.gz.tbi ${map.patient}_${map.sample}.candidateSV.vcf.gz.tbi

        mv manta/results/variants/diploidSV.vcf.gz ${map.patient}_${map.sample}.diploidSV.vcf.gz
        mv manta/results/variants/diploidSV.vcf.gz.tbi ${map.patient}_${map.sample}.diploidSV.vcf.gz.tbi

        mv manta/results/variants/somaticSV.vcf.gz ${map.patient}_${map.sample}.somaticSV.vcf.gz
        mv manta/results/variants/somaticSV.vcf.gz.tbi ${map.patient}_${map.sample}.somaticSV.vcf.gz.tbi
        """
    else
        """
        configManta.py \
        --normalBam ${map.normal} \
        --tumorBam ${map.tumor} \
        --referenceFasta ${params.ref} \
        --runDir manta 

        python2 manta/runWorkflow.py -j 8
        
        mv manta/results/variants/candidateSmallIndels.vcf.gz ${map.patient}_${map.sample}.candidateSmallIndels.vcf.gz
        mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi ${map.patient}_${map.sample}.candidateSmallIndels.vcf.gz.tbi

        mv manta/results/variants/candidateSV.vcf.gz ${map.patient}_${map.sample}.candidateSV.vcf.gz
        mv manta/results/variants/candidateSV.vcf.gz.tbi ${map.patient}_${map.sample}.candidateSV.vcf.gz.tbi

        mv manta/results/variants/diploidSV.vcf.gz ${map.patient}_${map.sample}.diploidSV.vcf.gz
        mv manta/results/variants/diploidSV.vcf.gz.tbi ${map.patient}_${map.sample}.diploidSV.vcf.gz.tbi
        
        mv manta/results/variants/somaticSV.vcf.gz ${map.patient}_${map.sample}.somaticSV.vcf.gz
        mv manta/results/variants/somaticSV.vcf.gz.tbi ${map.patient}_${map.sample}.somaticSV.vcf.gz.tbi

        """
}
