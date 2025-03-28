nextflow.enable.dsl=2

process MergeVcfs {
  conda "${params.java_env}"
  scratch true
  conda "${params.java_env}"
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(map), path(vcf)

  output:
  tuple val(map), path("*vcf"), emit: MUTECT2_vcf
  tuple val(map), path("*.idx"), emit: MUTECT2_idx

  script:
  def cmd = "java -jar ${params.database_path}/EVC_nextflow/picard/build/libs/picard.jar MergeVcfs"

  for( int i=0; i<20; i++ ) {
    cmd += " I= "
    cmd += vcf[i]
    }
  cmd += " O= ${map.patient}_${map.tumor_meta.sample}_mutect2_unfiltered.vcf"

  cmd

}
