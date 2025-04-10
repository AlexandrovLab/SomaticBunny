nextflow.enable.dsl=2

process Delly_SVcalling {
    scratch true
    conda "${params.delly_env}"
    label 'process_medium'
    publishDir("${params.Delly_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map), path("*.bcf"), path("*.bcf.csi"), emit: SVcalling_bcf

    script:
    """
    delly call -x ${params.database_dir}/Delly/human.hg38.excl.tsv  \
        -o ${map.patient}_${map.sample}.bcf \
        -g ${params.ref} \
        ${map.tumor} \
        ${map.normal} \
        -q 20 \
        -s 15 \
        -z 5
    """

}
