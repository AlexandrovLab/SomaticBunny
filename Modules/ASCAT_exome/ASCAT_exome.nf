nextflow.enable.dsl=2

process ASCAT_exome {
    conda "${params.ascat_env}"
    scratch true
    label 'process_low'
    publishDir("${params.ascat_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(patient), val(sample), val(gender), path(tumor_logr), path(tumor_baf), path(normal_logr), path(normal_baf)

    output:
    tuple val(patient), val(sample), path("Before_correction_*png"), emit: before_correction
    tuple val(patient), val(sample), path("After_correction_*png"), emit: after_correction
    tuple val(patient), val(sample), path("*ASPCF.png"), emit: segments_plot
    tuple val(patient), val(sample), path("*.txt"), emit: ASACT_txt_files
    tuple val(patient), val(sample), path("*.ASCATprofile.png"), emit: profile
    tuple val(patient), val(sample), path("*.sunrise.png"), emit: sunrise
    tuple val(patient), val(sample), path("purity_ploidy_*.txt"), emit: purityploidy

    script:
    """
    #!/usr/bin/env Rscript
    
    library(ASCAT)    
    
    ascat.bc = ascat.loadData(Tumor_LogR_file = "${tumor_logr}", Tumor_BAF_file = "${tumor_baf}", Germline_LogR_file = "${normal_logr}", Germline_BAF_file = "${normal_baf}", gender = "${gender}", genomeVersion = "hg38")
    
    ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
    ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "${params.database_dir}/ASCAT/WES/hg38/GC_Correction/GC_G1000_hg38.txt", replictimingfile = "${params.database_dir}/ASCAT/WES/hg38/RT_Correction/RT_G1000_hg38.txt")
    ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
    ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA)
    ascat.plotSegmentedData(ascat.bc)
    ascat.output = ascat.runAscat(ascat.bc, gamma = 1, write_segments = T)
    write.table(ascat.output[["segments"]], file = "${patient}_${sample}.segments.txt", sep = "\t", quote = F, row.names = F)
    QC = ascat.metrics(ascat.bc, ascat.output)
    write.table(QC, file = "${patient}_${sample}.metrics.txt", sep = "\t", quote = F, row.names = F)
    
    cnvs = ascat.output[["segments"]][2:6]
    write.table(cnvs, file = "${patient}_${sample}.cnvs.txt", sep = "\t", quote = F, row.names = F, col.names = T)
    
    # Extract and save purity and ploidy information
    purity <- ascat.output\$purity
    ploidy <- ascat.output\$ploidy
    
    # Create a data frame to store the information
    output_df <- data.frame(sample = "${patient}_${sample}", purity = purity, ploidy = ploidy)
    
    # Write the data frame to a text file
    write.table(output_df, file = "purity_ploidy_${patient}_${sample}.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = !file.exists("purity_ploidy_${patient}_${sample}.txt"))
    
    """
}
