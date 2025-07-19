nextflow.enable.dsl=2

workflow SAVE_CSV_Mutect2 {
    take:
        Mutect2_out
        outdir
        Mutect2_dir
    main:
        Mutect2_out.collectFile(keepHeader: true, storeDir: "${outdir}/csv") { patient -> 
            real_vcf = "${Mutect2_dir}/${patient}_mutect2_filtered.vcf"
            ["Mutect2.csv", "patient,vcf\n${patient},${real_vcf}\n"]
        }.set{Mutect2_file}
    
    emit:
        Mutect2_file

}