nextflow.enable.dsl=2

process POST {
    conda "${params.dkfz_env}"
    scratch true
    label 'process_low'  
    publishDir("${params.post_dir}/${patient_sample}/snvs", mode: 'copy', pattern: 'snvs_filtered/2outof4/*_PASSed.vcf')
    publishDir("${params.post_dir}/${patient_sample}/snvs", mode: 'copy', pattern: 'snvs_filtered/2outof4/*_snv_final_annotated.vcf')
    publishDir("${params.post_dir}/${patient_sample}/indels", mode: 'copy', pattern: 'indels_filtered/2outof3/*_PASSed.vcf')
    publishDir("${params.post_dir}/${patient_sample}/indels", mode: 'copy', pattern: 'indels_filtered/2outof3/*_indel_final_annotated.vcf')
    
    input:
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai)

    output:
    path("snvs_filtered/2outof4/${patient_sample}_PASSed.vcf"), emit: snv_passed
    path("snvs_filtered/2outof4/${patient_sample}_snv_final_annotated.vcf"), emit: snv_annotated, optional: true
    path("indels_filtered/2outof3/${patient_sample}_PASSed.vcf"), emit: indel_passed
    path("indels_filtered/2outof3/${patient_sample}_indel_final_annotated.vcf"), emit: indel_annotated, optional: true

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

    # ${params.database_path}/EVC_nextflow/Databases/filtering.sh

    bash /tscc/lustre/restricted/alexandrov-ddn/users/tiy002/projects/SMURFS/filtering.sh \
        ${patient_sample} \
        ${mutect2_vcf} \
        ${muse_vcf} \
        ${strelka_snv} \
        ${strelka_indel} \
        ${sage_vcf} \
        ${tumor_bam} \
        ${params.database_path} \
        ${params.ref}
    
    # Debug: Check output
    echo "=== Output Files ==="
    echo "Main directory:"
    ls -lh *_snv.vcf 2>/dev/null || echo "No SNV files in main directory"
    
    echo ""
    echo "SNVs filtered directory:"
    ls -lh snvs_filtered/
    wc -l snvs_filtered/*_snv.vcf 2>/dev/null || echo "No SNV files in snvs_filtered/"
    
    echo ""
    echo "SNV caller subdirectories:"
    ls -lh snvs_filtered/mutect_snvs/ 2>/dev/null || echo "No mutect_snvs directory"
    ls -lh snvs_filtered/muse_snvs/ 2>/dev/null || echo "No muse_snvs directory"
    ls -lh snvs_filtered/strelka_snvs/ 2>/dev/null || echo "No strelka_snvs directory"
    ls -lh snvs_filtered/sage_snvs/ 2>/dev/null || echo "No sage_snvs directory"
    
    echo ""
    echo "2outof4 directory:"
    ls -lh snvs_filtered/2outof4/
    
    echo ""
    echo "Final SNV files:"
    wc -l snvs_filtered/2outof4/${patient_sample}_snv_final_annotated.vcf || echo "File not found"
    wc -l snvs_filtered/2outof4/${patient_sample}_PASSed.vcf || echo "File not found"
    
    echo ""
    echo "Indels (for comparison):"
    wc -l indels_filtered/2outof3/${patient_sample}_indel_final_annotated.vcf || echo "File not found"
    wc -l indels_filtered/2outof3/${patient_sample}_PASSed.vcf || echo "File not found"
    """
}
