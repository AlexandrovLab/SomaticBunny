process RECALIBRATE_MergeBam {
    scratch true
    label 'RECALIBRATE'
    conda "${params.samtools_env}"
    errorStrategy = 'retry'
    maxRetries 3
        
    input:
    val(map)

    output:
    tuple val(map.patient), val(map.meta), path("*_temp_merged.bam"), emit: MergeBam_input

    script:
    def temp_merged = "${map.patient}_${map.meta.sample}_${map.status}_temp_merged.bam"
    
    """
    # Validate we have exactly 20 chunks
    if [ ${map.bam.size()} -ne 20 ]; then
        echo "ERROR: Expected 20 chunks, got ${map.bam.size()}"
        exit 1
    fi
    
    # Check all BAM files exist
    for bam in ${map.bam.join(' ')}; do
        if [ ! -f "\$bam" ]; then
            echo "ERROR: BAM file \$bam does not exist"
            exit 1
        fi
    done
    
    # Merge all 20 chunks (this removes duplicates from your current problematic merge)
    samtools merge -o ${temp_merged} ${map.bam.join(' ')}
    """
}
