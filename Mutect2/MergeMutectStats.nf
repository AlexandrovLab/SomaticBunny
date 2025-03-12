nextflow.enable.dsl=2

process MergeMutectStats {
  conda "${params.java_env}"
  scratch true
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
  maxRetries 3

  input:
  tuple val(patient), path(stats)

  output:
  tuple val(patient), path("*.stats"), emit: MUTECT2_stats

  script:
  def cmd = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk MergeMutectStats"

  for( int i=0; i<20; i++ ) {
    cmd += " -stats "
    cmd += stats[i]
    }
  cmd += " -O ${patient}_merged.stats"

  cmd

}
