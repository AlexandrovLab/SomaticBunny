nextflow.enable.dsl=2

process ASCAT_allelecount {
    conda "${params.alleleCounter_env}"
    scratch true
    label 'process_low'
    publishDir("${params.ascat_dir}", mode: 'copy')

    input:
    val(map)
    each chr

    output:
    tuple val(map), val(chr),
          path("${map.patient}_${map.sample}_tumor_${chr}.txt"),
          path("${map.patient}_${map.sample}_normal_${chr}.txt"), emit: allelecount

    script:
    """
    alleleCounter \
    -l "${params.ascat_loci_prefix_wes}${chr}.txt" \
    -b "${map.tumor}" \
    -o "${map.patient}_${map.sample}_tumor_${chr}.txt" -m 20

    alleleCounter \
    -l "${params.ascat_loci_prefix_wes}${chr}.txt" \
    -b "${map.normal}" \
    -o "${map.patient}_${map.sample}_normal_${chr}.txt" -m 20

    """
}
