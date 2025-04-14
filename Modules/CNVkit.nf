nextflow.enable.dsl=2

process CNVkit {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_high'
    publishDir("${params.cnvkit_dir}", mode: 'copy')
    // errorStrategy 'retry'
    // maxRetries 3

    input:
    val(map)
    path(reference_cnn)

    output:
    path("*.bed"), emit: CNVkit_bed
    path("*.cn*"), emit: CNVkit_cn_files
    path("*.pdf"), emit: CNVkit_pdf
    path("*.png"), emit: CNVkit_png

    script:
    if (params.type == "exome")
        """
        cnvkit.py batch \
        ${map.tumor} \
        -r ${reference_cnn} \
        -p 16 \
        --scatter --diagram

        cnvkit.py call \
        "${map.patient}_${map.sample}_tumor_recal.cns" \
        -o ${map.patient}_${map.sample}_calls.cns
        """
    else
        """
        cnvkit.py batch \
        ${tumor} \
        -r ${reference_cnn} \
        --method wgs \
        -p 16 \
        --scatter --diagram

        cnvkit.py call \
        "${map.patient}_${map.sample}_tumor_recal.cns" \
        -o ${map.patient}_${map.sample}_calls.cns
        """
}
