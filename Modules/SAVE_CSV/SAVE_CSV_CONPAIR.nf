nextflow.enable.dsl=2

workflow SAVE_CSV_CONPAIR {
    take:
        CONPAIR_out
        outdir
        CONPAIR_dir
    main:
        CONPAIR_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { patient, path -> 
            real_txt = "${CONPAIR_dir}/${patient}_info.txt"
            ["CONPAIR.csv", "patient,txt\n${patient},${real_txt}\n"]
        }.set{CONPAIR_file}
    
    emit:
        CONPAIR_file

}