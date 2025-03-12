nextflow.enable.dsl=2

workflow SAVE_CSV_MuSE2 {
    take:
        MuSE2_out
        outdir
        MuSE2_dir
    main:
        MuSE2_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { map -> 
            patient = map.patient
            real_vcf = "${MuSE2_dir}/${map.patient}.vcf"

            ["MuSE2.csv", "patient,vcf\n${patient},${real_vcf}\n"]
        }.set{MuSE2_file}
    
    emit:
        MuSE2_file

}