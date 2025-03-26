nextflow.enable.dsl=2

workflow SAVE_CSV_STRELKA {
    take:
        STRELKA_out
        outdir
        STRELKA_dir
    main:
        STRELKA_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { map -> 
            patient = map.patient
            real_vcf = "${STRELKA_dir}/${map.patient}.somatic_snvs.vcf.gz"

            ["STRELKA.csv", "patient,vcf\n${patient},${real_vcf}\n"]
        }.set{STRELKA_file}
    
    emit:
        STRELKA_file

}