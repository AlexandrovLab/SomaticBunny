nextflow.enable.dsl=2

process RECALIBRATE_MergeBam {
    scratch true
    label 'RECALIBRATE'
    conda "${params.samtools_env}"
    publishDir("${params.recal_dir}", mode: 'copy')
    errorStrategy = 'retry'
    maxRetries 3
        
    input:
    val(map)

    output:
    tuple val(map.patient), val(map.status), val(map.meta), path("*_temp_merged.bam"), emit: MergeBam_input


    script:
    def temp_merged = "${map.patient}_${map.meta.sample}_${map.status}_temp_merged.bam"
    
    // Build the merge command to temporary file
    def cmd = "samtools merge -o ${temp_merged}"

    for( int i=0; i<20; i++ ) {
        cmd += " ${map.bam[i]} "
    }

    cmd
}


