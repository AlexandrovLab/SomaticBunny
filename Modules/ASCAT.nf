nextflow.enable.dsl=2

process ASCAT {
    conda "${params.ascat_env}"
    scratch true
    label 'process_high'
    publishDir("${params.ascat_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.3.0--h41abebc_0' :
        'biocontainers/cancerit-allelecount:4.3.0--h41abebc_0' }"
        
    input:
    val(map)

    output:
    tuple val(map.meta), path("*alleleFrequencies*.txt"),      emit: allelefreq
    tuple val(map.meta), path("*png"),                             emit: pngs
    tuple val(map.meta), path("*BAF.txt"),                         emit: baf
    tuple val(map.meta), path("*LogR.txt"),                        emit: logr
    tuple val(map.meta), path("*raw.txt"),                        emit: cnvs
    tuple val(map.meta), path("*segments.txt"),                    emit: segments
    tuple val(map.meta), path("*cnvs.txt"),                        emit: cnv
    tuple val(map.meta), path("*metrics.txt"),                     emit: metrics
    tuple val(map.meta), path("purity_ploidy*.txt"),                emit: purityploidy


    script:
    """
    #!/usr/bin/env Rscript

    install.packages("BiocManager", repos="https://cloud.r-project.org")
    install.packages("ASCAT", repos="https://bioconductor.org/packages/release/bioc")

    library(ASCAT)

    # PrepareHTS: Extracting logR and BAF from HTS data (bam files)
    ascat.prepareHTS(
    tumourseqfile = "${map.tumor}",
    normalseqfile = "${map.normal}",
    tumourname = "${map.patient}_tumor",
    normalname = "${map.patient}_normal",
    allelecounter_exe = "alleleCounter",
    alleles.prefix = "${params.database_dir}/ASCAT/WGS/hg38/Alleles/G1000_alleles_hg38_chr",
    loci.prefix = "${params.database_dir}/ASCAT/WGS/hg38/Loci/G1000_loci_hg38_chr",
    gender = "${map.gender}",
    nthreads = "16",
    genomeVersion = "hg38")

    # Running ASCAT
    # For HTS data (WGS, WES and targeted sequencing), gamma must be set to 1 in ascat.runASCAT

    ascat.bc = ascat.loadData(Tumor_LogR_file = paste0("${map.patient}","_tumor_tumourLogR.txt"), 
                            Tumor_BAF_file = paste0("${map.patient}","_tumor_tumourBAF.txt"), 
                            Germline_LogR_file = paste0("${map.patient}","_tumor_normalLogR.txt"), 
                            Germline_BAF_file = paste0("${map.patient}","_tumor_normalBAF.txt"), 
                            gender = "${map.gender}", 
                            genomeVersion = "hg38")


    ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
    ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "${params.database_dir}/ASCAT/WGS/hg38/GC_Correction/GC_G1000_hg38.txt", replictimingfile = "${params.database_dir}/ASCAT/WGS/hg38/RT_Correction/RT_G1000_hg38.txt")
    ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
    ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA) # penalty=25 for targeted sequencing data
    ascat.plotSegmentedData(ascat.bc)
    ascat.output = ascat.runAscat(ascat.bc, gamma = 1, write_segments = T) # gamma=1 for HTS data
    write.table(ascat.output[["segments"]], file=paste0("${map.patient}", ".segments.txt"), sep="\t", quote=F, row.names=F)
    QC = ascat.metrics(ascat.bc,ascat.output)
    write.table(QC, file=paste0("${map.patient}", ".metrics.txt"), sep="\t", quote=F, row.names=F)

    # Save ASCAT objects
    save(ascat.bc, ascat.output, file = paste0("${map.patient}", "_ASCAT_objects.Rdata"))

    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("${map.patient}",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    # Extract and save purity and ploidy information
    purity <- ascat.output\$purity
    ploidy <- ascat.output\$ploidy

    # Create a data frame to store the information
    output_df <- data.frame(
        sample = "${map.patient}",
        purity = purity,
        ploidy = ploidy
    )

    # Write the data frame to a text file
    write.table(output_df, 
                file = "purity_ploidy_${map.patient}.txt", 
                append = TRUE, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = !file.exists("purity_ploidy_${map.patient}.txt"))

    """
}
