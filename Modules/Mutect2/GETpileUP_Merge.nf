nextflow.enable.dsl=2

process GETpileUP_Merge {
  scratch true
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
  maxRetries 3

  input:
  val(map)


  output:
  tuple val(map.patient), val(map.status), path("*table"), emit: CalculateContamination_input

  script:
  def cmd = "cat ${map.mix[0][0]}"

  for( int i=1; i<20; i++ ) {
    cmd += " <(tail -n +3 ${map.mix[i][0]}) "
    }
  cmd += " > ${map.patient}_getpileupsummaries_${map.status}.table"

  cmd

}