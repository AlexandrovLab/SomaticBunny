nextflow.enable.dsl=2

process BQSR_INTERVALS {
    conda "${params.gatk_env}"
    label 'process_low'
    publishDir("${params.intervals_dir}", mode: 'copy')

    input:
    path(packaged_dir)
    path(ref_dict)

    output:
    path "interval_list", emit: interval_dir

    script:
    """
    mkdir -p interval_list
    cp ${packaged_dir}/*-scattered.interval_list interval_list/
    python3 ${projectDir}/bin/make_contig_groups.py ${ref_dict} 20 interval_list/bqsr
    """
}
