nextflow.enable.dsl=2

process CNVkit {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_high'
    publishDir("${params.cnvkit_dir}", mode: 'copy')
    // errorStrategy 'retry'
    // maxRetries 3

    input:
    tuple val(patient), val(type), path(normal), path(tumor), path(cnn)

    output:
    path("*.bed"), emit: CNVkit_bed
    path("*.cn*"), emit: CNVkit_cn_files
    path("*.pdf"), emit: CNVkit_pdf
    path("*.png"), emit: CNVkit_png

    script:

    if (type == "exome")
        """
        cnvkit.py batch \
        ${tumor} \
        -r ${cnn} \
        -p 16 \
        --scatter --diagram

        cnvkit.py call \
        "${patient}_tumor_recal.cns" \
        -o ${patient}_calls.cns
        """
    else
        """
        cnvkit.py batch \
        ${tumor} \
        -r ${cnn} \
        --method wgs \
        -p 16 \
        --scatter --diagram

        cnvkit.py call \
        "${patient}_tumor_recal.cns" \
        -o ${patient}_calls.cns
        """
}
