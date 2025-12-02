process CLEANUP_BWA {
    label 'process_low'
    executor 'local' 
    
    input:
    val(meta)
    
    output:
    val(true)
    
    script:
    """
    # Clean up BWA_MEM work directories for this sample
    WORK_DIR="${workflow.workDir}"
    if [ -d "\$WORK_DIR" ]; then
        echo "Searching for BWA_MEM work directories for ${meta.patient}_${meta.sample}_${meta.status}..."
        
        # Find and delete BWA_MEM work directories
        find "\$WORK_DIR" -type f -name ".command.sh" -exec grep -l "bwa-mem2 mem.*${meta.patient}_${meta.sample}_${meta.status}" {} \\; 2>/dev/null | while read cmdfile; do
            workdir=\$(dirname "\$cmdfile")
            if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                echo "  Removing BWA_MEM work directory: \$workdir"
                rm -rf "\$workdir"
            fi
        done
        
        # Also remove raw BAM files from MKDUP work directories
        echo "Removing staged raw BAM from MKDUP directories for ${meta.patient}_${meta.sample}_${meta.status}..."
        find "\$WORK_DIR" -type f -name ".command.sh" -exec grep -l "picard MarkDuplicates.*${meta.patient}_${meta.sample}_${meta.status}" {} \\; 2>/dev/null | while read cmdfile; do
            workdir=\$(dirname "\$cmdfile")
            if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                raw_bam="\$workdir/${meta.patient}_${meta.sample}_${meta.status}_raw.bam"
                if [ -f "\$raw_bam" ]; then
                    echo "  Removing staged raw BAM: \$raw_bam"
                    rm -f "\$raw_bam"
                fi
            fi
        done
        
        echo "Cleanup completed for ${meta.patient}_${meta.sample}_${meta.status}"
    fi
    """
}
