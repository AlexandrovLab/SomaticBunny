nextflow.enable.dsl=2

workflow SAVE_CSV_MuSE {
    take:
        MuSE_out
        outdir
        MuSE_dir
    main:
        MuSE_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { map -> 
            patient = map.patient
            real_vcf = "${MuSE_dir}/${map.patient}.vcf"

            ["MuSE.csv", "patient,vcf\n${patient},${real_vcf}\n"]
        }.set{MuSE_file}
    
    emit:
        MuSE_file

}
