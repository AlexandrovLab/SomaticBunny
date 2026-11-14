nextflow.enable.dsl=2

process MKDUP {
    conda "${params.mkdup_env}"
    scratch true
    label 'MKDUP'
    if (params.publish.tokenize(',').contains('markdup_bam')) {
        publishDir("${params.mkdup_dir}", mode: 'copy')
    }
    maxRetries 3

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta.patient), val(meta), path("*bam"), path("*bai"), emit: pair_mutect
    tuple val(meta), path("*bam"), path("*bai"), emit: mkdup_bam
    val(meta), emit: cleanup_trigger

    script:
    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    // Spill files go to the task directory, which with `scratch true` is the node-local scratch space
    // the scheduler gives the job (TSCC: /scratch/$USER/job_$SLURM_JOB_ID), not to shared storage.
    // --mkdup_temp_dir overrides this for nodes whose local disk is too small for them.
    def tmp = params.mkdup_temp_dir ?: '$PWD/tmp'
    """
    mkdir -p ${tmp}
    picard MarkDuplicates \\
        -Xmx${xmx}M ${gc} -Djava.io.tmpdir=${tmp} \\
        TMP_DIR=${tmp} \\
        ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES=4000 MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.mkdup_optical_pixel_distance} \\
        VALIDATION_STRINGENCY=STRICT \\
        I=${bam} \\
        O=${meta.patient}_${meta.sample}_${meta.status}_mkdp.bam \\
        M=${meta.patient}_${meta.sample}_${meta.status}_markDuplicates_Matrix.txt
    """
}
