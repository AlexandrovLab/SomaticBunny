nextflow.enable.dsl=2

process POSTEVC {
    scratch true
    label 'process_low'
    publishDir("${params.postevc_dir}", mode: 'copy', saveAs: { 
            fn ->
             { "${postevc_map[0]}" }})

    input:
    val(postevc_map)

    output:
    path("*")

    script:
    """
    bash ${params.database_path}/EVC_nextflow/Databases/filtering.sh ${postevc_map[0]} ${postevc_map[1]} ${postevc_map[2]} ${postevc_map[4]} ${postevc_map[3]} ${postevc_map[5]} ${postevc_map[6]}
    """
}
