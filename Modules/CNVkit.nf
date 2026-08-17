nextflow.enable.dsl=2

process CNVkit {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_high'
    publishDir("${params.cnvkit_dir}", mode: 'copy')

    input:
    val map
    path reference_cnn

    output:
    path("*.bed"), emit: CNVkit_bed
    path("*.cn*"), emit: CNVkit_cn_files
    path("*.pdf"), emit: CNVkit_pdf
    path("*.png"), emit: CNVkit_png

    script:
    def tumor_prefix = map.tumor
        .toString()
        .tokenize("/")
        .last()
        .replaceFirst(/\.bam$/, "")

    def output_prefix = "${map.patient}_${map.sample}"
    def method_option = params.type == "exome" ? "" : "--method wgs"

    """
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH

    cnvkit.py batch \
        "${map.tumor}" \
        -r "${reference_cnn}" \
        ${method_option} \
        -p ${task.cpus} \
        --scatter \
        --diagram

    test -s "${tumor_prefix}.cns"

    cnvkit.py call \
        "${tumor_prefix}.cns" \
        -o "${output_prefix}_calls.cns"

    cnvkit.py export bed \
        "${output_prefix}_calls.cns" \
        -o "${output_prefix}_calls.bed"
    """
}