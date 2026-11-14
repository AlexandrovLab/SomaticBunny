
process RENAME_BAM_HEADER {
    tag "${meta.patient}_${meta.sample}"
    label 'process_low'
    conda "${params.samtools_env}"

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}.bam"), path("${meta.patient}_${meta.sample}_${meta.status}.bam.bai"), emit: renamed_bam

    script:
    def new_name   = "${meta.patient}_${meta.sample}_${meta.status}"
    def output_bam = "${new_name}.bam"
    """
    n_rg=\$(samtools view -H ${bam} | grep -c '^@RG' || true)

    if [ "\$n_rg" -gt 0 ]; then
        # Keep the existing read groups (ID, LB, PU, PL per lane/library); only SM becomes the pipeline sample name
        samtools view -H --no-PG ${bam} \\
        | awk -F'\\t' -v sm='${new_name}' 'BEGIN{OFS="\\t"}
            /^@RG/ { has=0; for (i=2;i<=NF;i++) if (\$i ~ /^SM:/) { \$i="SM:" sm; has=1 }
                     if (!has) \$0 = \$0 OFS "SM:" sm }
            { print }' > new_header.sam
        samtools reheader new_header.sam ${bam} > ${output_bam}
    else
        # No read group in the input: add a single one
        samtools addreplacerg -@ ${task.cpus} -m overwrite_all \\
            -r ID:${new_name} -r SM:${new_name} -r LB:${new_name} -r PL:ILLUMINA -r PU:${new_name} \\
            -o ${output_bam} ${bam}
    fi

    samtools index ${output_bam}
    """
}