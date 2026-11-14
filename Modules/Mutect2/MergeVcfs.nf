nextflow.enable.dsl = 2

process MergeVcfs {
    scratch true
    conda "${params.picard_merge_env}"
    label 'process_low'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(map), path(vcf)

    output:
    tuple val(map), path("*.vcf"), emit: MUTECT2_vcf
    tuple val(map), path("*.idx"), emit: MUTECT2_idx

    script:
    if (!(vcf instanceof List)) {
        throw new IllegalStateException(
            "MergeVcfs expected a list of 20 VCF shards, but received: ${vcf}"
        )
    }

    if (vcf.size() != 20) {
        throw new IllegalStateException(
            "MergeVcfs expected exactly 20 VCF shards for " +
            "${map.patient}_${map.tumor_meta.sample}, " +
            "but received ${vcf.size()}"
        )
    }

    def xmx = task.memory ? (task.memory.mega * 0.8).intValue() : 4096
    def gc  = "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${Math.max(1, task.cpus.intdiv(4))}"
    def cmd = "picard MergeVcfs -Xmx${xmx}M ${gc}"

    for (int i = 0; i < 20; i++) {
        cmd += " I=${vcf[i]}"
    }

    cmd += " O=${map.patient}_${map.tumor_meta.sample}_mutect2_unfiltered.vcf"

    cmd
}
