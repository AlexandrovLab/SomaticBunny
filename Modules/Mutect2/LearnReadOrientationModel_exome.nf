nextflow.enable.dsl=2

process LearnReadOrientationModel_exome {
  scratch true
  conda "${params.java_env}"
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(map), path(f1r2_files)

  output:
  tuple val(map), path("*read-orientation-model.tar.gz"), emit: MUTECT2_read_orientation

  script:
  """
  ${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk LearnReadOrientationModel -I ${f1r2_files} -O ${map.patient}_${map.tumor_meta.sample}_read-orientation-model.tar.gz
  """

}