nextflow.enable.dsl=2

process MergeMutectStats {
  scratch true
  conda "${params.java_env}"
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(map), path(stats)

  output:
  tuple val(map), path("*.stats"), emit: MUTECT2_stats

  script:
  def cmd = "${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk MergeMutectStats"

  for( int i=0; i<20; i++ ) {
    cmd += " -stats "
    cmd += stats[i]
    }
  cmd += " -O ${map.patient}_${map.tumor_meta.sample}_merged.stats"

  cmd

}