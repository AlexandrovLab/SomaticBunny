process RENAME_BAM_HEADER {
    tag "${meta.patient}_${meta.sample}"
    label 'RENAME_BAM_HEADER'
    conda "${params.mkdup_env}"
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}.bam"), path("${meta.patient}_${meta.sample}_${meta.status}.bam.bai"), emit: renamed_bam
    
    script:
    def new_name = "${meta.patient}_${meta.sample}_${meta.status}"
    def output_bam = "${new_name}.bam"
    """
    picard -Xmx${task.memory.toGiga()}g AddOrReplaceReadGroups \\
        I=${bam} \\
        O=${output_bam} \\
        RGID=${new_name} \\
        RGSM=${meta.patient}_${meta.sample}_${meta.status} \\
        RGPL=ILLUMINA \\
        RGLB=${meta.sample} \\
        RGPU=${new_name} \\
        SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        USE_JDK_DEFLATER=true \\
        USE_JDK_INFLATER=true
    
    mv ${new_name}.bai ${output_bam}.bai
    """
}
