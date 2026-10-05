nextflow.enable.dsl=2

process POST_scenario1 {
    conda "${params.POST_env}"
    scratch true
    label 'process_low'  
    publishDir("${params.post_dir}/${patient_sample}/scenario1", mode: 'copy')
    input:
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai)

    output:
    path("*.vcf")
    script:
    """
    set -e  # Exit on error
    
    # Create symlink so pysam can find the index file
    ln -sf ${tumor_bai} ${tumor_bam}.bai

    # Debug: Check conda environment
    echo "=== Conda Environment ==="
    which python
    python --version
    
    # Debug: Check input files
    echo "=== Input Files ==="
    ls -lh ${mutect2_vcf}
    ls -lh ${muse_vcf}
    ls -lh ${strelka_snv}
    ls -lh ${sage_vcf}
    ls -lh ${tumor_bam}

    bash ${projectDir}/scenario1.sh \
        ${patient_sample} \
        ${mutect2_vcf} \
        ${muse_vcf} \
        ${strelka_snv} \
        ${strelka_indel} \
        ${sage_vcf} \
        ${tumor_bam} \
        ${params.database_path} \
        ${params.ref}
    """

}
