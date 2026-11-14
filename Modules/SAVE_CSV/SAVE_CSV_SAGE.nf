nextflow.enable.dsl=2

workflow SAVE_CSV_SAGE {
    take:
        SAGE_out
        outdir
        SAGE_dir
    main:
        SAGE_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { map -> 
            patient = map.patient
            real_vcf = "${SAGE_dir}/${map.patient}.sage.vcf.gz"

            ["SAGE.csv", "patient,vcf\n${patient},${real_vcf}\n"]
        }.set{SAGE_file}
    
    emit:
        SAGE_file

}