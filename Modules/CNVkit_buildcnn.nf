nextflow.enable.dsl=2

process CNVkit_buildcnn {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_high'
    publishDir("${params.cnvkit_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path normal_files

    output:
    path("reference.cnn"), emit: CNVkit_ref_cnn
    path("*coverage.cnn"), emit: CNVkit_buildcnns

    script:
    def input_files = normal_files instanceof List ?
        normal_files : [normal_files]

    def normal_bams = input_files
        .findAll { it.toString().endsWith(".bam") }

    if (normal_bams.isEmpty()) {
        error "CNVkit_buildcnn received no normal BAM files"
    }

    def bam_files = normal_bams
        .collect { "\"${it}\"" }
        .join(" ")

    if (params.type == "exome") {
        """
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH

        cnvkit.py batch \
            -n ${bam_files} \
            --targets "${params.mosdepth_bed}" \
            --fasta "${params.ref}" \
            --output-reference reference.cnn \
            -p ${task.cpus}
        """
    } else {
        """
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH

        cnvkit.py batch \
            -n ${bam_files} \
            --method wgs \
            --fasta "${params.ref}" \
            --output-reference reference.cnn \
            -p ${task.cpus}
        """
    }
}