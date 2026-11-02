nextflow.enable.dsl=2

process ALLELECOUNTER {
    conda "${params.alleleCounter_env}"
    label 'process_low'
    publishDir(
        path: { "${params.VAF_dir}" },
        mode: 'copy'
    )
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.3.0--h41abebc_0' :
        'biocontainers/cancerit-allelecount:4.3.0--h41abebc_0' }"

    input:
    tuple val(map), val(scenario)

    output:
    tuple val(map.patient_sample), path("*alleleCounter.txt"),                          emit: allele_counts
    tuple val(map.patient_sample), path("${map.patient_sample}_${scenario}_VAF.tsv"),   emit: vaf_table

    script:
    def sample = "${map.patient_sample}_${scenario}"
    def ac_header = "#CHR\\tPOS\\tCount_A\\tCount_C\\tCount_G\\tCount_T\\tGood_depth"
    """
    # 1. unique SNV loci of the consensus VCF (2 columns, sorted): alleleCounter rejects '#' lines and
    #    would count a duplicated locus twice
    python3 ${params.vaf_script} loci --vcf ${map.final_vcf} --out ${sample}_snv_loci.tsv

    # 2. allele counts per SNV locus (alleleCounter defaults: base quality >= 20, mapping quality >= 35)
    if [ -s ${sample}_snv_loci.tsv ]; then
        alleleCounter -r ${params.ref} -l ${sample}_snv_loci.tsv -b ${map.tumor}  -o ${sample}_tumor_alleleCounter.txt
        alleleCounter -r ${params.ref} -l ${sample}_snv_loci.tsv -b ${map.normal} -o ${sample}_normal_alleleCounter.txt
    else
        printf '${ac_header}\\n' > ${sample}_tumor_alleleCounter.txt
        printf '${ac_header}\\n' > ${sample}_normal_alleleCounter.txt
    fi

    # 3. one table: SNV counts joined on (CHROM, POS) with one-to-one validation, indels counted from the BAMs
    #    by exact event after left-normalisation (same MAPQ / base-quality thresholds as alleleCounter)
    python3 ${params.vaf_script} table \\
        --vcf ${map.final_vcf} \\
        --tumor-counts ${sample}_tumor_alleleCounter.txt \\
        --normal-counts ${sample}_normal_alleleCounter.txt \\
        --tumor-bam ${map.tumor} --normal-bam ${map.normal} \\
        --ref ${params.ref} \\
        --sample ${sample} \\
        --out ${sample}_VAF.tsv
    """
}
