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
    val(map)

    output:
    tuple val(map.patient_sample), path("*alleleCounter.txt"),              emit: allele_counts
    tuple val(map.patient_sample), path("${map.patient_sample}_VAF.vcf"),   emit: vaf_vcf
    tuple val(map.patient_sample), path("${map.patient_sample}_vaf.vcf"),   emit: vaf_final

    script:
    patient_sample = map.patient_sample

    """

    # ── alleleCounter for SNV VAF ─────────────────────────────────────────────
    alleleCounter -r ${params.ref} -l ${map.final_vcf} -b ${map.tumor}  -o ${map.patient_sample}_tumor_alleleCounter.txt
    alleleCounter -r ${params.ref} -l ${map.final_vcf} -b ${map.normal} -o ${map.patient_sample}_normal_alleleCounter.txt

    sample=${patient_sample}
    normal=\${sample}_normal_alleleCounter.txt
    tumor=\${sample}_tumor_alleleCounter.txt
    vcf=${map.final_vcf}

    # ── SNV VAF via alleleCounter ─────────────────────────────────────────────
    grep -v "^#" \$vcf | sort -V -k1,1 -k2,2 > a
    cat \$tumor  | tail -n +2 | sort -V -k1,1 -k2,2 | cut -f3- > b
    cat \$normal | tail -n +2 | sort -V -k1,1 -k2,2 | cut -f3- > c

    paste -d "\t" a b c > d

    cat d | awk '{OFS="\t"
    if(\$14==0){\$20=\$21=\$22="no_AC_annotation"}
    else if(length(\$4)!=1){\$20=\$21=\$22="INDEL_SKIP"}
    else if(length(\$5)!=1){\$20=\$21=\$22="INDEL_SKIP"}
    else if(\$4=="A"){if(\$5=="C"){\$20=\$10":"\$11":"\$14;\$21=\$15":"\$16":"\$19;\$22=\$11/\$14}
                else if(\$5=="G"){\$20=\$10":"\$12":"\$14;\$21=\$15":"\$17":"\$19;\$22=\$12/\$14}
                else if(\$5=="T"){\$20=\$10":"\$13":"\$14;\$21=\$15":"\$18":"\$19;\$22=\$13/\$14}}
    else if(\$4=="C"){if(\$5=="A"){\$20=\$11":"\$10":"\$14;\$21=\$16":"\$15":"\$19;\$22=\$10/\$14}
                else if(\$5=="G"){\$20=\$11":"\$12":"\$14;\$21=\$16":"\$17":"\$19;\$22=\$12/\$14}
                else if(\$5=="T"){\$20=\$11":"\$13":"\$14;\$21=\$16":"\$18":"\$19;\$22=\$13/\$14}}
    else if(\$4=="G"){if(\$5=="A"){\$20=\$12":"\$10":"\$14;\$21=\$17":"\$15":"\$19;\$22=\$10/\$14}
                else if(\$5=="C"){\$20=\$12":"\$11":"\$14;\$21=\$17":"\$16":"\$19;\$22=\$11/\$14}
                else if(\$5=="T"){\$20=\$12":"\$13":"\$14;\$21=\$17":"\$18":"\$19;\$22=\$13/\$14}}
    else if(\$4=="T"){if(\$5=="A"){\$20=\$13":"\$10":"\$14;\$21=\$18":"\$15":"\$19;\$22=\$10/\$14}
                else if(\$5=="C"){\$20=\$13":"\$11":"\$14;\$21=\$18":"\$16":"\$19;\$22=\$11/\$14}
                else if(\$5=="G"){\$20=\$13":"\$12":"\$14;\$21=\$18":"\$17":"\$19;\$22=\$12/\$14}}
                tumor=\$20;normal=\$21;vaf=\$22;\$20="RD:AD:TD";\$21=tumor;\$22=normal;\$23=vaf
                print \$0}' | cut -f1-9,20- > e

    sed -e "1i\\#CHR\tLOC\tSAMPLE\tREF\tALT\tCALLERS\tFILTER\tDOT\tINFO\tFORMAT\tTUMOR\tNORMAL\tVAF" e > \${sample}_VAF.vcf

    # SNV final: drop INDEL_SKIP and no_AC_annotation rows
    awk 'NR == 1 {print "sample", \$0; next;}{print FILENAME, \$13;}' \${sample}_VAF.vcf > \${sample}_1.vcf
    cat \${sample}_1.vcf | tail -n +2 \
        | sed 's/_VAF.vcf//g' \
        | sed '/INDEL_SKIP/d' \
        | sed '/no_AC_annotation/d' \
        | awk -v OFS="\t" '\$1=\$1' > \${sample}_snv_vaf.vcf

    # ── INDEL VAF via pysam ───────────────────────────────────────────────────
    grep -v "^#" \$vcf | awk 'length(\$4) != 1 || length(\$5) != 1' > indels_raw.vcf || true

    if [ -s indels_raw.vcf ]; then
    python3 ${params.indel_vaf_script} \
        indels_raw.vcf \
        ${map.tumor} \
        ${map.normal} \
        \${sample} \
        \${sample}_indel_vaf.vcf
    else
        printf "sample\t#CHR\tLOC\tSAMPLE\tREF\tALT\tCALLERS\tFILTER\tDOT\tINFO\tFORMAT\tTUMOR\tNORMAL\tVAF\n" \
            > \${sample}_indel_vaf.vcf
    fi

    grep "^#" \${sample}_VAF.vcf > \${sample}_VAF_updated.vcf
    grep -v "^#" \${sample}_VAF.vcf | grep -v "INDEL_SKIP" >> \${sample}_VAF_updated.vcf
    tail -n +2 \${sample}_indel_vaf.vcf | cut -f2- >> \${sample}_VAF_updated.vcf
    grep "^#" \${sample}_VAF_updated.vcf > \${sample}_VAF.vcf
    grep -v "^#" \${sample}_VAF_updated.vcf | sort -k1,1V -k2,2n >> \${sample}_VAF.vcf
    rm -f \${sample}_VAF_updated.vcf

    # ── Merge SNV + INDEL into final _vaf.vcf ────────────────────────────────
    cat \${sample}_snv_vaf.vcf > \${sample}_vaf.vcf
    tail -n +2 \${sample}_indel_vaf.vcf >> \${sample}_vaf.vcf

    # ── Clean up ──────────────────────────────────────────────────────────────
    rm -f a b c d e \${sample}_1.vcf indels_raw.vcf
    """
}
