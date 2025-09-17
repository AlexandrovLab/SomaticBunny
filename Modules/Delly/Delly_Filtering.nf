nextflow.enable.dsl=2

process Delly_Filtering {
    scratch true
    conda "${params.delly_env}"
    label 'process_medium'
    publishDir("${params.Delly_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(map), path(pre_bcf), path(csi)
    path normal_files

    output:
    tuple val(map), path("*geno.bcf"), path("*.csi"), emit: Delly_geno_bcf

    script:
    def base_name = pre_bcf.simpleName
    base_name = base_name.replace(".pre", "")
    """
    echo "Processing file: ${pre_bcf}"
    echo "Base name extracted: ${base_name}"
    
    delly call -g ${params.ref} \
        -v ${pre_bcf}  \
        -o ${base_name}_geno.bcf  \
        -x ${params.delly_excl}  \
        ${map.tumor} \
        ${normal_files.findAll { it.toString().endsWith('.bam') }.join(' ')} \
        -q 20 \
        -s 15 \
        -z 5

    # Create index for the output BCF file
    bcftools index ${base_name}_geno.bcf
    """
}
