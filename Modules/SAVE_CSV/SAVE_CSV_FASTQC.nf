nextflow.enable.dsl=2

workflow SAVE_CSV_FASTQC {
    take:
        FASTQC_out
        outdir
        FASTQC_dir
    main:
        FASTQC_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { meta, html -> 
            patient = meta.patient
            status = meta.status
            fastq_1 = meta.fastq_1
            fastq_2 = meta.fastq_2
            type = meta.type
            real_html = "${FASTQC_dir}/${html}"
            ["FASTQC.csv", "patient,status,fastq_1,fastq_2,type,html\n${patient},${status},${fastq_1},${fastq_2},${type},${real_html}\n"]
        }.set{FASTQC_file}
    
    emit:
        FASTQC_file

}