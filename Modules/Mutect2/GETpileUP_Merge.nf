nextflow.enable.dsl=2

process GETpileUP_Merge {
  scratch true
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy 'retry'
  maxRetries 3

  input:
  val(map)

  output:
  tuple val(map.patient), val(map.meta), path("*table"), emit: CalculateContamination_input

  script:
  def cmd = "cat ${map.mix[0][0]}"

  for( int i=1; i<20; i++ ) {
    cmd += " <(tail -n +3 ${map.mix[i][0]}) "
    }
  cmd += " > ${map.patient}_getpileupsummaries_${map.meta.sample}_${map.meta.status}.table"

  cmd

}
