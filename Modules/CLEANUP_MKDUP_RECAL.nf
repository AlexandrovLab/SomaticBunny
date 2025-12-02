process CLEANUP_MKDUP_RECAL {
    label 'process_low'
    executor 'local'
    
    input:
    val(meta)
    
    output:
    val(true)
    
    script:
    """
    # Give time for any file operations to complete
    sleep 5
    
    WORK_DIR="${workflow.workDir}"
    SAMPLE_ID="${meta.patient}_${meta.sample}_${meta.status}"
    
    echo "Starting cleanup for \${SAMPLE_ID}..."
    
    if [ -d "\$WORK_DIR" ]; then
        # 1. Clean up mkdup BAM files
        echo "Cleaning up mkdup BAM files..."
        find "\$WORK_DIR" -type f -name "\${SAMPLE_ID}_mkdp.bam*" 2>/dev/null | while read file; do
            workdir=\$(dirname "\$file")
            if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                echo "  Removing: \$file"
                rm -f "\$file"
            fi
        done
        
        # 2. Clean up intermediate recalibrated chunk BAM files
        echo "Cleaning up intermediate recalibrated chunk BAM files..."
        find "\$WORK_DIR" -type f -name "\${SAMPLE_ID}_recalibrated_[0-9]*.bam*" 2>/dev/null | while read file; do
            workdir=\$(dirname "\$file")
            if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                echo "  Removing: \$file"
                rm -f "\$file"
            fi
        done
        
        # 3. Clean up temporary merged BAM files (from RECALIBRATE_MergeBam)
        echo "Cleaning up temporary merged BAM files..."
        find "\$WORK_DIR" -type f -name "\${SAMPLE_ID}_temp_merged.bam*" 2>/dev/null | while read file; do
            workdir=\$(dirname "\$file")
            if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                echo "  Removing: \$file"
                rm -f "\$file"
            fi
        done
        
        # 4. Clean up MKDUP work directories
        echo "Cleaning up MKDUP work directories..."
        find "\$WORK_DIR" -type f -name ".command.sh" 2>/dev/null | while read cmdfile; do
            if grep -q "picard MarkDuplicates.*\${SAMPLE_ID}_mkdp.bam" "\$cmdfile" 2>/dev/null; then
                workdir=\$(dirname "\$cmdfile")
                if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                    if [ ! -f "\$workdir/.lock" ]; then
                        echo "  Removing MKDUP work directory: \$workdir"
                        rm -rf "\$workdir"
                    fi
                fi
            fi
        done
        
        # 5. Clean up RECALIBRATE_BQSR work directories
        echo "Cleaning up RECALIBRATE_BQSR work directories..."
        find "\$WORK_DIR" -type f -name ".command.sh" 2>/dev/null | while read cmdfile; do
            if grep -q "gatk ApplyBQSR.*\${SAMPLE_ID}_recalibrated_[0-9]" "\$cmdfile" 2>/dev/null; then
                workdir=\$(dirname "\$cmdfile")
                if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                    if [ ! -f "\$workdir/.lock" ]; then
                        echo "  Removing RECALIBRATE_BQSR work directory: \$workdir"
                        rm -rf "\$workdir"
                    fi
                fi
            fi
        done
        
        # 6. Clean up RECALIBRATE_MergeBam work directories
        echo "Cleaning up RECALIBRATE_MergeBam work directories..."
        find "\$WORK_DIR" -type f -name ".command.sh" 2>/dev/null | while read cmdfile; do
            if grep -q "samtools merge.*\${SAMPLE_ID}_temp_merged.bam" "\$cmdfile" 2>/dev/null; then
                workdir=\$(dirname "\$cmdfile")
                if [ -f "\$workdir/.exitcode" ] && [ "\$(cat "\$workdir/.exitcode")" = "0" ]; then
                    if [ ! -f "\$workdir/.lock" ]; then
                        echo "  Removing RECALIBRATE_MergeBam work directory: \$workdir"
                        rm -rf "\$workdir"
                    fi
                fi
            fi
        done
        
        echo "Cleanup completed for \${SAMPLE_ID}"
    fi
    """
}