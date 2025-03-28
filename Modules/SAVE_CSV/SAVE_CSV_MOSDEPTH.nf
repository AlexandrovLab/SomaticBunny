nextflow.enable.dsl=2

workflow SAVE_CSV_MOSDEPTH {
    take:
        MOSDEPTH_out
        outdir
        MOSDEPTH_dir
    main:
        MOSDEPTH_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { patient, path -> 
            real_txt = "${MOSDEPTH_dir}/${patient}_tumor.mosdepth.summary.txt"
            ["MOSDEPTH.csv", "patient,vcf\n${patient},${real_txt}\n"]
        }.set{MOSDEPTH_file}
    
    emit:
        MOSDEPTH_file

}