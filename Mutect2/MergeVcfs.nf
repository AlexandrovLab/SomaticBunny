nextflow.enable.dsl=2

process MergeVcfs {
  conda "${params.java_env}"
  scratch true
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
  maxRetries 3

  input:
  tuple val(patient), path(vcf)

  output:
  tuple val(patient), path("*vcf"), emit: MUTECT2_vcf
  tuple val(patient), path("*.idx"), emit: MUTECT2_idx

  script:
  def cmd = "java -jar /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/picard/build/libs/picard.jar MergeVcfs"

  for( int i=0; i<20; i++ ) {
    cmd += " I= "
    cmd += vcf[i]
    }
  cmd += " O= ${patient}_mutect2_unfiltered.vcf"

  cmd

}
