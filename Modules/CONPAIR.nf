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
    tuple val(map.patient), path("*txt"), emit: CONPAIR_out

    script:
    """
    # Use a task-local temporary directory.
    mkdir -p "\$PWD/tmp_conpair"
    export TMPDIR="\$PWD/tmp_conpair"

    # Force the pinned Conda Java instead of a Java installation inherited
    # from the user's TSCC environment.
    export JAVA_HOME="\$CONDA_PREFIX"
    export PATH="\$CONDA_PREFIX/bin:\$PATH"
    hash -r

    # Java settings for GATK3.
    unset _JAVA_OPTIONS
    export JAVA_TOOL_OPTIONS="-Xms4g -Xmx16g -Djava.io.tmpdir=\$TMPDIR"

    # External Conpair v0.2 source directory.
    export CONPAIR_DIR="${params.conpair}"
    export PYTHONPATH="\${PYTHONPATH:-}:${params.conpair}/modules"

    # Locate the GATK3 JAR installed by the pinned Conda environment.
    GATK_JAR=\$(find "\$CONDA_PREFIX" \
        -type f \
        -name 'GenomeAnalysisTK.jar' \
        -print \
        -quit)

    if [[ -z "\$GATK_JAR" || ! -f "\$GATK_JAR" ]]; then
        echo "ERROR: GenomeAnalysisTK.jar was not found in \$CONDA_PREFIX" >&2
        exit 1
    fi

    echo "Python: \$(command -v python)"
    echo "Java:   \$(command -v java)"
    echo "GATK:   \$GATK_JAR"

    \$CONDA_PREFIX/bin/python ${params.conpair}/scripts/run_gatk_pileup_for_sample.py \
    -B ${map.normal} \
    -O ${map.patient}_${map.tumor_meta.sample}_normal.pileup \
    -D ${params.conpair} \
    -G "\$GATK_JAR" \
    --reference ${params.ref} \
    --markers ${params.conpair_marker};
    
    \$CONDA_PREFIX/bin/python ${params.conpair}/scripts/run_gatk_pileup_for_sample.py \
    -B ${map.tumor} \
    -O ${map.patient}_${map.tumor_meta.sample}_tumor.pileup \
    -D ${params.conpair} \
    -G "\$GATK_JAR" \
    --reference ${params.ref} \
    --markers ${params.conpair_marker};

    awk -F" " '\$5!=""' ${map.patient}_${map.tumor_meta.sample}_normal.pileup > ${map.patient}_${map.tumor_meta.sample}_normal.cleanpileup;
    awk -F" " '\$5!=""' ${map.patient}_${map.tumor_meta.sample}_tumor.pileup > ${map.patient}_${map.tumor_meta.sample}_tumor.cleanpileup;

    \$CONDA_PREFIX/bin/python ${params.conpair}/scripts/estimate_tumor_normal_contamination.py \
    -T ${map.patient}_${map.tumor_meta.sample}_tumor.cleanpileup \
    -N ${map.patient}_${map.tumor_meta.sample}_normal.cleanpileup \
    --outfile ${map.patient}_${map.tumor_meta.sample}_contamination.txt \
    --markers ${params.conpair_marker_txt};
    
    \$CONDA_PREFIX/bin/python ${params.conpair}/scripts/verify_concordance.py \
    -T ${map.patient}_${map.tumor_meta.sample}_tumor.cleanpileup \
    -N ${map.patient}_${map.tumor_meta.sample}_normal.cleanpileup \
    --normal_homozygous_markers_only \
    --min_cov 10 \
    --outfile ${map.patient}_${map.tumor_meta.sample}_concordance.txt \
    --markers ${params.conpair_marker_txt};
    
    echo "====Contamination Estimation" > ${map.patient}_${map.tumor_meta.sample}_info.txt;
    cat ${map.patient}_${map.tumor_meta.sample}_contamination.txt >> ${map.patient}_${map.tumor_meta.sample}_info.txt;
    echo "====Concordance" >> ${map.patient}_${map.tumor_meta.sample}_info.txt;
    cat ${map.patient}_${map.tumor_meta.sample}_concordance.txt >> ${map.patient}_${map.tumor_meta.sample}_info.txt;

    rm ${map.patient}_${map.tumor_meta.sample}_contamination.txt  ${map.patient}_${map.tumor_meta.sample}_concordance.txt
    """

}
