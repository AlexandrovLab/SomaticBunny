nextflow.enable.dsl=2

process RECALIBRATE_SortBam {
    scratch true
    label 'RECALIBRATE_SortBam'
    conda "${params.samtools_env}"
    if (params.publish.tokenize(',').contains('recal_bam')) {
        publishDir("${params.recal_dir}", mode: 'copy')
    }
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(mkdup_bam), path(mkdup_bai)

    output:
    tuple val(patient), val(meta), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam"), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam.bai"), emit: pair_recal
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam"), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam.bai"), emit: recal_out
    val(meta), emit: cleanup_trigger

    script:
    def sorted_bam = "${meta.patient}_${meta.sample}_${meta.status}_recal.bam"
    """
    samtools sort -@ ${task.cpus} ${bam} -o ${sorted_bam}
    samtools index ${sorted_bam}

    # The ApplyBQSR scatter/gather must neither lose nor duplicate reads: per-contig record counts (mapped + unmapped,
    # '*' = unplaced unmapped) of the recalibrated BAM must equal those of the MarkDuplicates BAM. Index only, instant.
    samtools idxstats ${mkdup_bam} | awk -F'\\t' '{print \$1"\\t"\$3+\$4}' | sort > expected.counts
    samtools idxstats ${sorted_bam} | awk -F'\\t' '{print \$1"\\t"\$3+\$4}' | sort > observed.counts
    if ! diff expected.counts observed.counts > counts.diff; then
        echo "ERROR: read counts per contig differ between ${mkdup_bam} (<) and ${sorted_bam} (>): reads were lost or duplicated in the ApplyBQSR scatter/gather"
        cat counts.diff
        exit 1
    fi
    echo "read count check passed: \$(awk -F'\\t' '{s+=\$2} END{print s}' observed.counts) records in both BAMs"
    """
}
