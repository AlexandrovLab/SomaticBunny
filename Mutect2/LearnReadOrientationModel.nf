nextflow.enable.dsl=2

process LearnReadOrientationModel {
  conda "${params.java_env}"
  scratch true
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
  maxRetries 3

  input:
  tuple val(patient), path(f1r2_files)


  output:
  tuple val(patient), path("*read-orientation-model.tar.gz"), emit: MUTECT2_read_orientation

  script:
  def cmd = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk LearnReadOrientationModel"

  for( int i=0; i<20; i++ ) {
    cmd += " -I "
    cmd += f1r2_files[i]
    }
  cmd += " -O ${patient}_read-orientation-model.tar.gz"

  cmd

}
