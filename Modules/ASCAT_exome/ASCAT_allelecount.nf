nextflow.enable.dsl=2

process ASCAT_allelecount {
    conda "${params.alleleCounter_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.ascat_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    val(map)
    each chr

    output:
    tuple val(map), val(chr), path("*_${chr}.txt"), emit: allelecount

    script:
    """
    alleleCounter \
    -l "${params.database_dir}/ASCAT/WES/hg38/Loci/G1000_loci_hg38_chr${chr}.txt" \
    -b "${map.tumor}" \
    -o "${map.patient}_${map.sample}_tumor_${chr}.txt" -m 20

    alleleCounter \
    -l "${params.database_dir}/ASCAT/WES/hg38/Loci/G1000_loci_hg38_chr${chr}.txt" \
    -b "${map.normal}" \
    -o "${map.patient}_${map.sample}_normal_${chr}.txt" -m 20

    """
}
