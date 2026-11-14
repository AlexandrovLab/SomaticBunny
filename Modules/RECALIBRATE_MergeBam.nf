nextflow.enable.dsl=2

/*
 * Gathers the 20 ApplyBQSR chunk BAMs of one sample into its recalibrated BAM.
 *
 * Each chunk holds whole contigs (bqsr/<n>-contigs.interval_list; chunk 20 also the unmapped
 * reads) and comes out of ApplyBQSR coordinate-sorted, so `samtools merge` produces the final
 * coordinate-sorted BAM directly and no re-sort is needed. -c and -p keep the one @RG and @PG
 * line the chunks share; without them samtools gives every chunk after the first a renamed copy.
 *
 * Also checks that the scatter/gather neither lost nor duplicated reads: per-contig record counts
 * of the merged BAM must equal those of the MarkDuplicates BAM it was computed from.
 */
process RECALIBRATE_MergeBam {
    scratch true
    label 'RECALIBRATE_MergeBam'
    conda "${params.samtools_env}"
    if (params.publish.tokenize(',').contains('recal_bam')) {
        publishDir("${params.recal_dir}", mode: 'copy')
    }

    input:
    tuple val(patient), val(meta), path(chunks), path(mkdup_bam), path(mkdup_bai)

    output:
    tuple val(patient), val(meta), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam"), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam.bai"), emit: pair_recal
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam"), path("${meta.patient}_${meta.sample}_${meta.status}_recal.bam.bai"), emit: recal_out
    val(meta), emit: cleanup_trigger

    script:
    def recal_bam = "${meta.patient}_${meta.sample}_${meta.status}_recal.bam"
    def n_chunks  = chunks instanceof List ? chunks.size() : 1
    if (n_chunks != 20) {
        throw new IllegalStateException(
            "RECALIBRATE_MergeBam expected 20 ApplyBQSR chunks for ${meta.patient}_${meta.sample}_${meta.status}, got ${n_chunks}")
    }
    """
    samtools merge -@ ${task.cpus} -c -p -o ${recal_bam} ${chunks.join(' ')}
    samtools index -@ ${task.cpus} ${recal_bam}

    # Per-contig record counts (mapped + unmapped, '*' = unplaced unmapped) must match the
    # MarkDuplicates BAM. Index only, instant.
    samtools idxstats ${mkdup_bam} | awk -F'\\t' '{print \$1"\\t"\$3+\$4}' | sort > expected.counts
    samtools idxstats ${recal_bam} | awk -F'\\t' '{print \$1"\\t"\$3+\$4}' | sort > observed.counts
    if ! diff expected.counts observed.counts > counts.diff; then
        echo "ERROR: read counts per contig differ between ${mkdup_bam} (<) and ${recal_bam} (>): reads were lost or duplicated in the ApplyBQSR scatter/gather"
        cat counts.diff
        exit 1
    fi
    echo "read count check passed: \$(awk -F'\\t' '{s+=\$2} END{print s}' observed.counts) records in both BAMs"
    """
}
