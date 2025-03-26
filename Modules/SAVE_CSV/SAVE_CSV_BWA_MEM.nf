nextflow.enable.dsl=2

workflow SAVE_CSV_BWA_MEM {
    take:
        BWA_MEM_out
        outdir
        bam_dir
    main:
        BWA_MEM_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { meta, bam -> 
            patient = meta.patient
            status = meta.status
            fastq_1 = meta.fastq_1
            fastq_2 = meta.fastq_2
            type = meta.type
            real_bam = "${bam_dir}/${meta.patient}_${meta.status}_raw.bam"
            ["BWA_MEM.csv", "patient,status,fastq_1,fastq_2,type,bam\n${patient},${status},${fastq_1},${fastq_2},${type},${real_bam}\n"]
        }.set{BWA_MEM_file}
    
    emit:
        BWA_MEM_file

}