nextflow.enable.dsl=2

process POST_rename {
    conda "${params.POST_rename_env}"
    scratch true
    label 'process_low' 
    input:
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai)

    output:
    tuple val(patient_sample), path("*_mutect.vcf"), path("*_muse.vcf"), path("*_strelka_snv.vcf"), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai), emit: renamed_postInput
    script:
    """
    suffix_1="mutect2"
    suffix_2="strelka"
    suffix_3="muse" 

    echo "INFO/DP DP_\$suffix_1" > rename_mutect.txt
    echo "INFO/DP DP_\$suffix_2" > rename_strelka.txt
    echo "INFO/SOMATIC SOMATIC_\$suffix_2" >> rename_strelka.txt
    echo "INFO/SOMATIC SOMATIC_\$suffix_3" > rename_muse.txt
    bcftools annotate --rename-annots rename_mutect.txt ${mutect2_vcf} -o ${patient_sample}rename_mutect.vcf
    bcftools annotate --rename-annots rename_strelka.txt ${strelka_snv} -o ${patient_sample}rename_strelka_snv.vcf 
    bcftools annotate --rename-annots rename_muse.txt ${muse_vcf} -o ${patient_sample}rename_muse.vcf 
    """
}
