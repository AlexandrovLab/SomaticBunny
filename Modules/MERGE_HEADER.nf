nextflow.enable.dsl=2

process MERGE_HEADER {
    conda "${params.POST_env}"
    scratch true
    label 'process_low'
    publishDir("${params.post_dir}/standard_vcf/${patient_sample}/*_standard_PASSed.vcf")

    input:
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai)
    tuple val(patient_sample), path(final_vcf), val(scenario)

    output:
    path("*_standard_PASSed.vcf"), emit: final_standard_vcf
    script:
    """
    set -e  # Exit on error
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

    bash ${projectDir}/rename_info.sh \
        ${mutect2_vcf} \
        ${strelka_snv} \
        ${muse_vcf}

    bash ${projectDir}/merge_header.sh \
        output_mutect.vcf \
        output_strelka.vcf \
        output_muse.vcf \
        ${sage_vcf} > merged_header.vcf

    cat merged_header.vcf ${final_vcf} > ${sample}_${scenario}_standard_PASSed.vcf
    """

}
