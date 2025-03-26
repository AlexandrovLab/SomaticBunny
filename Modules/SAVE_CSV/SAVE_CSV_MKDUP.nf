nextflow.enable.dsl=2

workflow SAVE_CSV_MKDUP {
    take:
        MKDUP_out
        outdir
        bam_dir
    main:
        MKDUP_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { meta, bam, bai -> 
            patient = meta.patient
            status = meta.status
            fastq_1 = meta.fastq_1
            fastq_2 = meta.fastq_2
            type = meta.type
            real_bam = "${bam_dir}/${meta.patient}_${meta.status}_mkdp.bam"
            real_bai = "${bam_dir}/${meta.patient}_${meta.status}_mkdp.bai"
            ["MKDUP.csv", "patient,status,fastq_1,fastq_2,type,bam,bai\n${patient},${status},${fastq_1},${fastq_2},${type},${real_bam},${real_bai}\n"]
        }.set{MKDUP_file}
    
    emit:
        MKDUP_file

}