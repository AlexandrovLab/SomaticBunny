nextflow.enable.dsl=2

process LearnReadOrientationModel {
  scratch true
  conda "${params.gatk_env}"
  label 'process_low'
  publishDir("${params.MUTECT2_dir}", mode: 'copy')
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(map), path(f1r2_files)

  output:
  tuple val(map), path("*read-orientation-model.tar.gz"), emit: MUTECT2_read_orientation

  script:
  def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
  def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
  def cmd = "${params.gatk} LearnReadOrientationModel --java-options \"-Xmx${xmx}M ${gc}\""

  for( int i=0; i<20; i++ ) {
    cmd += " -I "
    cmd += f1r2_files[i]
    }
  cmd += " -O ${map.patient}_${map.tumor_meta.sample}_read-orientation-model.tar.gz"

  cmd

}
