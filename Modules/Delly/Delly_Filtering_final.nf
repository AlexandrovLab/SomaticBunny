nextflow.enable.dsl=2

process Delly_Filtering_final {
    scratch true
    conda "${params.delly_env}"
    label 'process_medium'
    publishDir("${params.Delly_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(map), path(geno_bcf), path(csi), path(sample_tsv)

    output:
    tuple val(map), path("*somatic.bcf"), emit: Delly_somatic_bcf

    script:
    def base_name = geno_bcf.simpleName
    """
    delly filter -f somatic  \
	-o ${base_name}.somatic.bcf \
	-s ${sample_tsv}  \
	${geno_bcf}
    """

}
