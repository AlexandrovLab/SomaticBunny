nextflow.enable.dsl=2

process CONPAIR {
    conda "${params.conpair_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.conpair_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map.patient), path("*txt"), emit: coverage
    tuple val(map.patient), path("*txt"), emit: CONPAIR_out

    script:
    """
    export PATH=${params.jre}/bin:$PATH
    export CONPAIR_DIR=${params.conpair}
    export GATK_JAR=${params.database_dir}/GenomeAnalysisTK.jar
    export PYTHONPATH=\${PYTHONPATH:-}:${params.conpair}/modules/

    \$CONDA_PREFIX/bin/python2 ${params.conpair}/scripts/run_gatk_pileup_for_sample.py \
    -B ${map.normal} \
    -O ${map.patient}_normal.pileup \
    -D ${params.conpair} \
    -G ${params.database_dir}/GenomeAnalysisTK.jar \
    --reference ${params.ref} \
    --markers ${params.database_dir}/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed;
    
    \$CONDA_PREFIX/bin/python2 ${params.conpair}/scripts/run_gatk_pileup_for_sample.py \
    -B ${map.tumor} \
    -O ${map.patient}_tumor.pileup \
    -D ${params.conpair} \
    -G ${params.database_dir}/GenomeAnalysisTK.jar \
    --reference ${params.ref} \
    --markers ${params.database_dir}/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed;

    awk -F" " '\$5!=""' ${map.patient}_normal.pileup > ${map.patient}_normal.cleanpileup;
    awk -F" " '\$5!=""' ${map.patient}_tumor.pileup > ${map.patient}_tumor.cleanpileup;

    \$CONDA_PREFIX/bin/python2 ${params.conpair}/scripts/estimate_tumor_normal_contamination.py \
    -T ${map.patient}_tumor.cleanpileup \
    -N ${map.patient}_normal.cleanpileup \
    --outfile ${map.patient}_contamination.txt \
    --markers ${params.database_dir}/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt;
    
    \$CONDA_PREFIX/bin/python2 ${params.conpair}/scripts/verify_concordance.py \
    -T ${map.patient}_tumor.cleanpileup \
    -N ${map.patient}_normal.cleanpileup \
    --normal_homozygous_markers_only \
    --min_cov 10 \
    --outfile ${map.patient}_concordance.txt \
    --markers ${params.database_dir}/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt;
    
    echo "====Contamination Estimation" > ${map.patient}_info.txt;
    cat ${map.patient}_contamination.txt >> ${map.patient}_info.txt;
    echo "====Concordance" >> ${map.patient}_info.txt;
    cat ${map.patient}_concordance.txt >> ${map.patient}_info.txt;

    rm ${map.patient}_contamination.txt  ${map.patient}_concordance.txt
    """

}