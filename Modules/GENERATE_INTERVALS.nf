process GENERATE_INTERVALS {
    conda "${params.gatk_env}"
    scratch true
    label 'process_low'
    publishDir("${params.intervals_dir}", mode: 'copy')

    input:
    path(ref)
    path(ref_fai)
    path(ref_dict)

    output:
    path "interval_list/*-scattered.interval_list", emit: interval_files
    path "interval_list", emit: interval_dir

    script:
    """
    # Build chromosome list from .fai
    if grep -q "^chr" ${ref_fai}; then
        awk '\$1 ~ /^chr([0-9]+|X|Y|M)\$/' ${ref_fai} | cut -f1 > main_chromosomes.list
    else
        cut -f1 ${ref_fai} > main_chromosomes.list
    fi

    mkdir -p interval_list

    ${params.gatk} SplitIntervals \
        -R ${ref} \
        -L main_chromosomes.list \
        --scatter-count 20 \
        -O interval_list

    # Rename from 0-indexed to 1-indexed
    cd interval_list
    for f in *-scattered.interval_list; do
        idx=\$(echo "\${f}" | grep -oP '^\\d+')
        num=\$((10#\${idx} + 1))
        mv "\${f}" "\${num}-scattered.interval_list"
    done
    """
}
