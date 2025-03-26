nextflow.enable.dsl=2

process LearnReadOrientationModel {
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
  def cmd = "${params.database_path}/EVC_nextflow/gatk-4.6.0.0/gatk LearnReadOrientationModel"

  for( int i=0; i<20; i++ ) {
    cmd += " -I "
    cmd += f1r2_files[i]
    }
  cmd += " -O _read-orientation-model.tar.gz"

  cmd

}