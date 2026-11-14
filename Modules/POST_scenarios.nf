nextflow.enable.dsl=2

// Scenario 6 thresholds; override in nextflow.config or on the command line.
params.tlod        = 10
params.somatic_evs = 13
params.pon_veto    = false   // true: scenario 6 drops Mutect2 panel_of_normals sites

// Regions allowed into the somatic consensus (audit points 1 and 6).
// consensus_contigs: comma list or a file with one contig per line. callable_bed: BED of
// callable regions (e.g. the WES targets). With neither, human references use chr1-22, X, Y
// (or 1-22, X, Y); mm39 and RN7 use the primary chromosomes listed below. Everything else
// (viral/EBV, decoy, alt, unplaced, chrM) is kept out of the consensus.
params.consensus_contigs = null
params.callable_bed      = null

/*
 * Runs all six consensus scenarios in one task with consensus_scenarios.sh.
 *
 * Published to ${params.post_dir}/${patient_sample}/scenario_N/ (N = 1..6):
 *   scenario_N/${patient_sample}_final.vcf                  SNVs + indels, sorted, one header
 *   scenario_N/snvs/${patient_sample}_snv.vcf
 *   scenario_N/snvs/${patient_sample}_snv_annotated.tsv.gz
 *   scenario_N/indels/${patient_sample}_indel.vcf
 *   scenario_N/indels/${patient_sample}_indel_annotated.tsv.gz
 */
process POST_SCENARIOS {
    conda "${params.POST_env}"
    scratch true
    label 'process_low'

    publishDir "${params.post_dir}/${patient_sample}", mode: 'copy',
        saveAs: { fn ->
            def m = (fn =~ /^scenario([1-6])\/(.+)$/)
            if (!m.matches()) return null
            def n = m.group(1)
            def name = m.group(2)
            if (name == "${patient_sample}_final.vcf")                                   return "scenario_${n}/${name}"
            if (name.endsWith('_snv.vcf')   || name.endsWith('_snv_annotated.tsv.gz'))   return "scenario_${n}/snvs/${name}"
            if (name.endsWith('_indel.vcf') || name.endsWith('_indel_annotated.tsv.gz')) return "scenario_${n}/indels/${name}"
            return null
        }

    input:
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(strelka_indel), path(sage_vcf), path(tumor_bam), path(tumor_bai)

    output:
    // Scenario 1: all raw calls
    path("scenario1/${patient_sample}_snv.vcf"),   emit: snv_s1
    path("scenario1/${patient_sample}_indel.vcf"), emit: indel_s1
    tuple val(patient_sample), path("scenario1/${patient_sample}_final.vcf"), emit: final_vcf_s1
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(sage_vcf), path("scenario1/${patient_sample}_final.vcf"), emit: merge_input_s1

    // Scenario 2: PASS calls from any caller
    path("scenario2/${patient_sample}_snv.vcf"),   emit: snv_s2
    path("scenario2/${patient_sample}_indel.vcf"), emit: indel_s2
    tuple val(patient_sample), path("scenario2/${patient_sample}_final.vcf"), emit: final_vcf_s2
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(sage_vcf), path("scenario2/${patient_sample}_final.vcf"), emit: merge_input_s2

    // Scenario 3: >= 2 callers, any FILTER
    path("scenario3/${patient_sample}_snv.vcf"),   emit: snv_s3
    path("scenario3/${patient_sample}_indel.vcf"), emit: indel_s3
    tuple val(patient_sample), path("scenario3/${patient_sample}_final.vcf"), emit: final_vcf_s3
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(sage_vcf), path("scenario3/${patient_sample}_final.vcf"), emit: merge_input_s3

    // Scenario 4: >= 2 callers, at least one PASS
    path("scenario4/${patient_sample}_snv.vcf"),   emit: snv_s4
    path("scenario4/${patient_sample}_indel.vcf"), emit: indel_s4
    tuple val(patient_sample), path("scenario4/${patient_sample}_final.vcf"), emit: final_vcf_s4
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(sage_vcf), path("scenario4/${patient_sample}_final.vcf"), emit: merge_input_s4

    // Scenario 5: >= 2 PASS callers
    path("scenario5/${patient_sample}_snv.vcf"),   emit: snv_s5
    path("scenario5/${patient_sample}_indel.vcf"), emit: indel_s5
    tuple val(patient_sample), path("scenario5/${patient_sample}_final.vcf"), emit: final_vcf_s5
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(sage_vcf), path("scenario5/${patient_sample}_final.vcf"), emit: merge_input_s5

    // Scenario 6: S5 + bias filter + TLOD/SomaticEVS thresholds + PON
    path("scenario6/${patient_sample}_snv.vcf"),   emit: snv_s6
    path("scenario6/${patient_sample}_indel.vcf"), emit: indel_s6
    tuple val(patient_sample), path("scenario6/${patient_sample}_final.vcf"), emit: final_vcf_s6
    tuple val(patient_sample), path(mutect2_vcf), path(muse_vcf), path(strelka_snv), path(sage_vcf), path("scenario6/${patient_sample}_final.vcf"), emit: merge_input_s6

    // Annotated tables for all scenarios (also published per scenario)
    path("scenario*/${patient_sample}_*_annotated.tsv.gz"), emit: annotated

    script:
    def primary = [
        mm39: ((1..19).collect { "chr${it}" } + ['chrX', 'chrY']).join(','),
        RN7 : ((1..20).collect { "chr${it}" } + ['chrX', 'chrY']).join(',')
    ]
    def contigs = params.consensus_contigs ?: primary[params.genome]
    def region_args = [
        contigs             ? "-c ${contigs}"             : '',
        params.callable_bed ? "-b ${params.callable_bed}" : ''
    ].findAll().join(' ')
    """
    set -euo pipefail

    # Index for the bias filter (pysam) must sit next to the BAM as <bam>.bai.
    if [[ "${tumor_bai}" != "${tumor_bam}.bai" ]]; then
        ln -sfn "${tumor_bai}" "${tumor_bam}.bai"
    fi

    echo "=== Environment ==="
    which python bcftools
    python --version
    bcftools --version | head -1

    echo "=== Input files ==="
    ls -lhL ${mutect2_vcf} ${muse_vcf} ${strelka_snv} ${strelka_indel} ${sage_vcf} ${tumor_bam} ${tumor_bam}.bai

    bash ${projectDir}/consensus_scenarios.sh \\
        -s all \\
        ${region_args} \\
        -t ${params.tlod} \\
        -e ${params.somatic_evs} \\
        ${params.pon_veto ? '-p' : ''} \\
        ${patient_sample} \\
        ${mutect2_vcf} \\
        ${muse_vcf} \\
        ${strelka_snv} \\
        ${strelka_indel} \\
        ${sage_vcf} \\
        ${tumor_bam} \\
        ${params.database_path} \\
        ${params.ref}

    # Per scenario: combine SNVs and indels into one sorted VCF with a single header.
    # Header lines from both files are merged (duplicates removed, order kept) so every
    # INFO key used by an SNV or an indel record is declared.
    echo "=== Output ==="
    for n in 1 2 3 4 5 6; do
        dir=scenario\${n}
        snv=\${dir}/${patient_sample}_snv.vcf
        indel=\${dir}/${patient_sample}_indel.vcf
        final=\${dir}/${patient_sample}_final.vcf
        {
            awk '/^##/ && !seen[\$0]++' "\${snv}" "\${indel}"
            awk '/^#CHROM/ { print; exit }' "\${snv}"
            awk '!/^#/' "\${snv}" "\${indel}" | sort -t\$'\\t' -k1,1V -k2,2n -k4,4 -k5,5
        } > "\${final}"

        # Fail the task if the combined file is not a readable VCF.
        n_final=\$(bcftools view -H "\${final}" | wc -l)
        n_snv=\$(awk '!/^#/' "\${snv}" | wc -l)
        n_indel=\$(awk '!/^#/' "\${indel}" | wc -l)
        if [[ \${n_final} -ne \$((n_snv + n_indel)) ]]; then
            echo "ERROR: \${final} has \${n_final} records, expected \$((n_snv + n_indel))" >&2
            exit 1
        fi
        echo "scenario \${n}: \${n_snv} SNVs, \${n_indel} indels -> \${final}"

        gzip -f \${dir}/${patient_sample}_snv_annotated.tsv \${dir}/${patient_sample}_indel_annotated.tsv
    done
    """
}
