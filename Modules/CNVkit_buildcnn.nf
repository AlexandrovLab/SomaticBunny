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
    def bam_files = normal_files.findAll { it.toString().endsWith('.bam') }.join(' ')
    
    if (params.type == "exome")
        """
        cnvkit.py batch \
        ${bam_files} \
        -n \
        --targets ${params.mosdepth_bed} \
        --fasta ${params.ref} \
        --output-reference reference.cnn \
        -p 8
        """
    else
        """
        cnvkit.py batch \
        ${bam_files} \
        -n \
        --method wgs \
        --fasta ${params.ref} \
        --output-reference reference.cnn \
        -p 8
        """
}

