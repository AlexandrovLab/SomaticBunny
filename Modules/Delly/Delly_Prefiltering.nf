nextflow.enable.dsl=2

process Delly_Prefiltering {
    scratch true
    conda "${params.delly_env}"
    label 'process_medium'
    publishDir("${params.Delly_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(map), path(bcf), path(csi), path(sample_tsv)

    output:
    tuple val(map), path("*pre.bcf"), path("*.csi"), emit: Delly_pre_bcf

    script:
    def base_name = bcf.simpleName
    """
    echo "Processing BCF file: ${bcf}"
    echo "Expected pattern based on map: ${map.patient}_${map.sample}.bcf"

    delly filter -f somatic  \
	-o ${base_name}.pre.bcf \
	-s ${sample_tsv}  \
	${bcf}
    """

}
