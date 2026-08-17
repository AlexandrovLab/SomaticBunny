// Author: Ting Yang, George Wu
// Lab: The Alexandrov Lab @ UCSD
// Date: 2026.6
// Version: 3.0

nextflow.enable.dsl=2

params.database_path = null
params.log_file = "$projectDir/pipeline.log" 
params.sample = "sample.csv"
params.genome = ""
params.tool = ""

// Initialize the log file with a header
new File(params.log_file).text = """
==============================================
Pipeline Run - ${new Date()}
==============================================

"""
// Validate database being set
if (!params.database_path) {
    error "ERROR: --database_path is required. Please provide the path to your SomaticBunny_database directory.\n" +
          "Usage: nextflow run main.nf --database_path /path/to/SomaticBunny_database"
}

params.genomes = [
    'GRCh38': [
        ref: "${params.database_path}/GRCh38_ref/GRCh38.d1.vd1.fa",
        bed: "${params.database_path}/GRCh38_ref/hg38_chr.bed.gz",
        conpair_marker: "${params.database_path}/Databases/GRCh38/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed",
        conpair_marker_txt: "${params.database_path}/Databases/GRCh38/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt",
        mosdepth_bed: "${params.database_path}/Databases/GRCh38/GRCh38_exome.bed",
        muse_vcf: "${params.database_path}/Databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        recal_knownsite1: "${params.database_path}/Databases/GRCh38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf",
        recal_knownsite2: "${params.database_path}/Databases/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz",
        recal_interval_wes: "${params.database_path}/Databases/GRCh38/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list",
        sage: [
            sage_ref_dir: "${params.database_path}/SAGE/v5_34/ref/38",
            sage_hotspots: "${params.database_path}/SAGE/v5_34/ref/38/variants/KnownHotspots.somatic.38.vcf.gz",
            sage_panel_bed: "${params.database_path}/SAGE/v5_34/ref/38/variants/ActionableCodingPanel.38.bed.gz",
            sage_high_confidence_bed: "${params.database_path}/SAGE/v5_34/ref/38/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz",
            sage_ref_genome_version: "38" 
        ],
        ASCAT: [
            ascat_ref_genome_version: "hg38",
            ascat_allele_prefix: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/Alleles/G1000_alleles_hg38_chr",
            ascat_loci_prefix: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/Loci/G1000_loci_hg38_chr",
            ascat_GCcontentfile: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/GC_Correction/GC_G1000_hg38.txt",
            ascat_replictimingfile: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/RT_Correction/RT_G1000_hg38.txt",
			ascat_allele_prefix_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/Alleles/G1000_alleles_hg38_chr",
			ascat_loci_prefix_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/Loci/G1000_loci_hg38_chr",
            ascat_GCcontentfile_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/GC_Correction/GC_G1000_hg38.txt",
            ascat_replictimingfile_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/RT_Correction/RT_G1000_hg38.txt"
        ],
        database_subdir: "GRCh38",
        mutect2_pon: "MuTect2.PON.5210.vcf.gz",
        mutect2_pon_wes: "MuTect2.PON.5210.vcf.gz",
        mutect2_germline: "af-only-gnomad.hg38_no_alt.vcf.gz",
        mutect2_interval_dir: "GRCh38_interval_list_20",
        delly_excl: "${params.database_path}/Databases/GRCh38/Delly/human.hg38.excl.tsv",
        tools: ["fastqc", "bwa_mem", "mkdup", "recalibrate", "sage", "strelka", "muse", "mutect2", "ascat", "delly", "cnvkit", "mosdepth", "conpair", "manta"]
    ],
    'GRCh37': [
        ref: "${params.database_path}/GRCh37_ref/GRCh37.fa", 
        bed: "${params.database_path}/GRCh37_ref/GRCh37_chr.bed.gz",
        conpair_marker: "${params.database_path}/Databases/GRCh37/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed",
        conpair_marker_txt: "${params.database_path}/Databases/GRCh37/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt",
        mosdepth_bed: "${params.database_path}/Databases/GRCh37/GRCh37_exome.bed",
        muse_vcf: "${params.database_path}/Databases/GRCh37/Homo_sapiens_assembly19.dbsnp.vcf",
        recal_knownsite1: "${params.database_path}/Databases/GRCh37/Homo_sapiens_assembly19.dbsnp.vcf.gz",
        recal_knownsite2: "${params.database_path}/Databases/GRCh37/Homo_sapiens_assembly19.known_indels.vcf.gz",
        recal_interval_wes: "${params.database_path}/Databases/GRCh37/intervals_Broad.human.exome.b37.interval_list",
        sage: [
            sage_ref_dir: "${params.database_path}/SAGE/v5_34/ref/37",
            sage_hotspots: "${params.database_path}/SAGE/v5_34/ref/37/dna/variants/KnownHotspots.somatic.37.vcf.gz",
            sage_panel_bed: "${params.database_path}/SAGE/v5_34/ref/37/dna/variants/ActionableCodingPanel.37.bed.gz",
            sage_high_confidence_bed: "${params.database_path}/SAGE/v5_34/ref/37/dna/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz",
            sage_ref_genome_version: "37" 
        ],
        ASCAT: [
            ascat_ref_genome_version: "hg19",
            ascat_allele_prefix: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/Alleles/G1000_alleles_hg19_chr",
            ascat_loci_prefix: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/Loci/G1000_loci_hg19_chr",
            ascat_GCcontentfile: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/GC_Correction/GC_G1000_hg19.txt",
            ascat_replictimingfile: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/RT_Correction/RT_G1000_hg19.txt",
            ascat_allele_prefix_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/Alleles/G1000_alleles_hg19_chr",
            ascat_loci_prefix_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/Loci/G1000_loci_hg19_chr",
            ascat_GCcontentfile_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/GC_Correction/GC_G1000_hg19.txt",
            ascat_replictimingfile_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/RT_Correction/RT_G1000_hg19.txt"
        ],
        database_subdir: "GRCh37",
        mutect2_pon: "Mutect2-WGS-panel-b37.vcf",
        mutect2_pon_wes: "Mutect2-exome-panel_b37.vcf",
        mutect2_germline: "af-only-gnomad.raw.sites.grch37.vcf.gz",
        mutect2_interval_dir: "GRCh37_interval_list_20",
        delly_excl: "${params.database_path}/Databases/GRCh37/Delly/human.hg19.excl.tsv",
        tools: ["fastqc", "bwa_mem", "mkdup", "recalibrate", "sage", "strelka", "muse", "mutect2", "ascat", "delly", "cnvkit", "mosdepth", "conpair", "manta"]
    ],
    'mm39': [
        ref: "${params.database_path}/mm39_ref/mm39.fa",
        bed: "${params.database_path}/mm39_ref/mm39_chr.bed.gz",
        database_subdir: "mm39",
        mosdepth_bed: "${params.database_path}/Databases/mm39/mm39_exome.bed",
        recal_interval_wes: "${params.database_path}/Databases/mm39/mm39_exome.interval_list",
        muse_vcf: "af_only_mgp_mm39_unique.vcf.gz",
        mutect2_pon: "PoN.mm39.vcf.gz",
        mutect2_pon_wes: "PoN.mm39.vcf.gz",
        mutect2_germline: "af-only-mgp.mm39.vcf.gz",
        mutect2_interval_dir: "mm39_interval_list_20",
        tools: ["fastqc", "bwa_mem", "mkdup", "strelka", "muse", "mutect2", "cnvkit", "mosdepth", "manta"]
    ],
    'RN7': [
        ref: "${params.database_path}/RN7_ref/rn7.fa",
        bed: "${params.database_path}/RN7_ref/rn7_chr.bed.gz",
        database_subdir: "RN7",
        mosdepth_bed: "${params.database_path}/Databases/RN7/rn7_exome.bed",
        recal_interval_wes: "${params.database_path}/Databases/RN7/rn7_exome.interval_list",
        muse_vcf: "",
        mutect2_pon: "",
        mutect2_pon_wes: "",
        mutect2_germline: "af-only_rn7.vcf.gz",
        mutect2_interval_dir: "RN7_interval_list_20",
        tools: ["fastqc", "bwa_mem", "mkdup", "strelka", "muse", "mutect2", "cnvkit", "mosdepth", "manta"]
    ]
]

// Custom ref/bed is REQUIRED when starting from markdup, recalibration, or variant_calling
if (params.first_step in ['markdup', 'recalibration', 'variant_calling']) {
    def missing = []
    if (!params.ref)      missing << '--ref'
    if (!params.ref_fai)  missing << '--ref_fai'
    if (!params.ref_dict) missing << '--ref_dict'
    if (!params.bed)      missing << '--bed'
    if (!params.bed_tbi)  missing << '--bed_tbi'
    if (missing) {
        error "When --first_step is '${params.first_step}', you must provide all custom reference files: " +
              "--ref, --ref_fai, --ref_dict, --bed, --bed_tbi. Missing: ${missing.join(', ')}"
    }
}

// If any custom ref/bed provided, validate they come as a complete set
if (params.ref || params.bed || params.bed_tbi || params.ref_fai || params.ref_dict) {
    def missing = []
    if (!params.ref)      missing << '--ref'
    if (!params.ref_fai)  missing << '--ref_fai'
    if (!params.ref_dict) missing << '--ref_dict'
    if (!params.bed)      missing << '--bed'
    if (!params.bed_tbi)  missing << '--bed_tbi'
    if (missing) {
        error "When providing custom reference files, all five must be specified: " +
              "--ref, --ref_fai, --ref_dict, --bed, --bed_tbi. Missing: ${missing.join(', ')}"
    }
    if (!file(params.ref).exists())      error "Custom ref file not found: ${params.ref}"
    if (!file(params.ref_fai).exists())  error "Custom ref_fai file not found: ${params.ref_fai}"
    if (!file(params.ref_dict).exists()) error "Custom ref_dict file not found: ${params.ref_dict}"
    if (!file(params.bed).exists())      error "Custom bed file not found: ${params.bed}"
    if (!file(params.bed_tbi).exists())  error "Custom bed_tbi file not found: ${params.bed_tbi}"

    log.info """
    ==============================================
    Using CUSTOM reference files:
      ref      : ${params.ref}
      ref_fai  : ${params.ref_fai}
      ref_dict : ${params.ref_dict}
      bed      : ${params.bed}
      bed_tbi  : ${params.bed_tbi}
    ==============================================
    """
}


// Validate genome selection
if (!params.genomes.containsKey(params.genome)) {
    error "Invalid genome: ${params.genome}. Available options are: ${params.genomes.keySet().join(', ')}"
}

if (params.custom_ref) {
    log.warn """
    WARNING: You are using a custom reference genome.
    Ensure your reference uses the same contig naming convention 
    (e.g., chr1 vs 1) as the pipeline's database files for genome '${params.genome}'.
    Mismatched contig names will cause tools like Mutect2, MuSE, BQSR, 
    SAGE, Conpair, and Delly to fail.
    """
}

// Check if user inputs custom ref or not
params.custom_ref = (params.ref && params.ref != params.genomes[params.genome].ref)


// Set genome-specific parameters
params.ref = params.ref ?: params.genomes[params.genome].ref
params.ref_fai = params.ref_fai ?: "${params.ref}.fai"
params.ref_dict = params.ref_dict ?: params.ref.replaceAll(/\.fa(sta)?$/, '.dict')
params.bed = params.bed ?: params.genomes[params.genome].bed
params.bed_tbi = params.bed_tbi ?: "${params.bed}.tbi"
params.conpair_marker = params.conpair_marker ?: params.genomes[params.genome].conpair_marker
params.conpair_marker_txt = params.conpair_marker_txt ?: params.genomes[params.genome].conpair_marker_txt
params.mosdepth_bed = params.mosdepth_bed ?: params.genomes[params.genome].mosdepth_bed
params.muse_vcf = params.muse_vcf ?: params.genomes[params.genome].muse_vcf
params.recal_knownsite1 = params.recal_knownsite1 ?: params.genomes[params.genome].recal_knownsite1
params.recal_knownsite2 = params.recal_knownsite2 ?: params.genomes[params.genome].recal_knownsite2
params.recal_interval_wes = params.recal_interval_wes ?: params.genomes[params.genome].recal_interval_wes
params.delly_excl = params.delly_excl ?: params.genomes[params.genome].delly_excl

// Set SAGE-specific parameters
params.SAGE_ref_dir = params.SAGE_ref_dir ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_ref_dir : null)
params.SAGE_hotspots = params.SAGE_hotspots ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_hotspots : null)
params.SAGE_panel_bed = params.SAGE_panel_bed ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_panel_bed : null)
params.SAGE_high_confidence_bed = params.SAGE_high_confidence_bed ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_high_confidence_bed : null)
params.SAGE_ref_genome_version = params.SAGE_ref_genome_version ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_ref_genome_version : null)

// Set ASCAT-specific parameters
params.ascat_ref_genome_version = params.ascat_ref_genome_version ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_ref_genome_version : null)
params.ascat_allele_prefix = params.ascat_allele_prefix ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_allele_prefix : null)
params.ascat_loci_prefix = params.ascat_loci_prefix ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_loci_prefix : null)
params.ascat_GCcontentfile = params.ascat_GCcontentfile ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_GCcontentfile : null)
params.ascat_replictimingfile = params.ascat_replictimingfile ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_replictimingfile : null)

params.ascat_allele_prefix_wes = params.ascat_allele_prefix_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_allele_prefix_wes : null)
params.ascat_loci_prefix_wes = params.ascat_loci_prefix_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_loci_prefix_wes : null)
params.ascat_GCcontentfile_wes = params.ascat_GCcontentfile_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_GCcontentfile_wes : null)
params.ascat_replictimingfile_wes = params.ascat_replictimingfile_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_replictimingfile_wes : null)


params.database_dir = "${params.database_path}/Databases"
params.genome_database_dir = "${params.database_dir}/${params.genomes[params.genome].database_subdir}"

// Genome-specific file paths for Mutect2
params.mutect2_pon = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_pon}"
params.mutect2_pon_wes = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_pon_wes}"
params.mutect2_interval_dir = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_interval_dir}"
params.mutect2_germline = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_germline}"

params.mutect2_targets = params.type == "exome" ? 
    "${params.mutect2_interval_dir}/${params.genomes[params.genome].mutect2_targets_exome}" : 
    "${params.mutect2_interval_dir}/${params.genomes[params.genome].mutect2_targets_wgs}"

// Define which tools are available for the selected genome
def available_tools = params.genomes[params.genome].tools

// Create a function to write to the log file
def writeToLog(logFile, message) {
    new File(logFile).append(message + "\n")
}

// Generate comprehensive configuration log
def configInfo = """
Pipeline Configuration:
==============================================
Run ID         : ${workflow.runName}
Command Line   : ${workflow.commandLine}
Profile        : ${workflow.profile}
Container      : ${workflow.container}
Run as User    : ${workflow.userName}
Launch Dir     : ${workflow.launchDir}
Work Dir       : ${workflow.workDir}
Project Dir    : ${projectDir}
==============================================
Analysis Parameters:
==============================================
Analysis Type  : ${params.type}
Genome         : ${params.genome}
==============================================
Reference Files:
==============================================
ref            : ${params.ref}
ref_fai        : ${params.ref_fai}
bed            : ${params.bed}
bed_tbi        : ${params.bed_tbi}
==============================================
Avaliable Tools:
==============================================
${available_tools.join(', ')}
==============================================
"""

params.bam_dir="$projectDir/RESULTS/BAM"
params.report_dir="$projectDir/RESULTS/REPORT"
params.manta_bed="${params.database_path}/GRCh38_ref/manta.bed.gz"
params.mkdup_temp_dir="$projectDir/mkdup_tmp"
params.FASTQC="${params.database_path}/FastQC"
params.FASTQC_dir="$projectDir/RESULTS/FASTQC"
params.mkdup_dir="$projectDir/RESULTS/MKDUP"
params.recal_dir="$projectDir/RESULTS/RECALIBRATE"
params.SAGE_jar="${params.database_path}/SAGE/sage_v3.3.jar"
params.SAGE_dir="$projectDir/RESULTS/SAGE"
params.strelka_dir="$projectDir/RESULTS/STRELKA"
params.muse_dir="$projectDir/RESULTS/MuSE"
params.MUTECT2_dir="$projectDir/RESULTS/Mutect2"
params.mosdepth_dir="$projectDir/RESULTS/mosdepth"
params.conpair_dir="$projectDir/RESULTS/Conpair"
params.cnvkit_dir="$projectDir/RESULTS/CNVkit"
params.Delly_dir="$projectDir/RESULTS/Delly"
params.ascat_dir="$projectDir/RESULTS/ASCAT"
params.conpair="${params.database_path}/Conpair-0.2"
params.manta_dir="$projectDir/RESULTS/MANTA"
params.post_dir="$projectDir/RESULTS/POST"
params.intervals_dir="$projectDir/RESULTS/custom_intervals"
params.VAF_dir="$projectDir/RESULTS/POST/VAF"

params.bwamem2_env = "$projectDir/yml/bwamem2.yml"
params.mkdup_env = "$projectDir/yml/mkdup.yml"
params.conpair_env = "$projectDir/yml/conpair.yml"
params.samtools_env = "$projectDir/yml/samtools.yml"
params.strelka_env = "$projectDir/yml/strelka_env.yml"
params.mosdepth_env = "$projectDir/yml/mosdepth_env.yml"
params.summary_env = "$projectDir/yml/py_summary.yml"
params.fastqc_env = "$projectDir/yml/fastqc_env.yml"
params.cnvkit_env = "$projectDir/yml/cnvkit.yml"
params.delly_env = "$projectDir/yml/delly.yml"
params.ascat_env = "$projectDir/yml/ascat.yml"
params.alleleCounter_env = "$projectDir/yml/alleleCounter.yml"
params.manta_env = "$projectDir/yml/manta.yml"
params.dkfz_env = "$projectDir/yml/dkfz.yml"
params.SAGE_java_env = "$projectDir/yml/sage_java.yml"
params.muse_env = "$projectDir/yml/muse.yml"
params.picard_merge_env = "$projectDir/yml/picard_merge.yml"
params.gatk_env = "${projectDir}/yml/gatk_runtime.yml"
params.gatk = "${params.database_path}/gatk-4.6.0.0/gatk"

params.tmp_dir = '${workflow.workDir}/ascat_tmp'

// Write configuration to log file
writeToLog(params.log_file, configInfo)

// Also display in terminal
log.info configInfo

include { FASTQC } from './Modules/FASTQC'
include { BWA_MEM } from './Modules/BWA_MEM'
include { CHECK_BAM_BWA } from './Modules/CHECK_BAM_BWA'
include { CLEANUP_BWA } from './Modules/CLEANUP_BWA'
include { COMBINE_REPORTS_BWA } from './Modules/COMBINE_REPORTS_BWA'
include { RECALIBRATE_BaseRecal } from './Modules/RECALIBRATE_BaseRecal.nf'
include { RECALIBRATE_BaseRecal_exome } from './Modules/RECALIBRATE_BaseRecal_exome.nf'
include { RECALIBRATE_MergeReport } from './Modules/RECALIBRATE_MergeReport.nf'
include { RECALIBRATE_BQSR } from './Modules/RECALIBRATE_BQSR.nf'
include { RECALIBRATE_BQSR_exome } from './Modules/RECALIBRATE_BQSR_exome.nf'
include { RECALIBRATE_MergeBam } from './Modules/RECALIBRATE_MergeBam.nf'
include { RECALIBRATE_SortBam } from './Modules/RECALIBRATE_SortBam.nf'
include { CHECK_BAM_RECAL } from './Modules/CHECK_BAM_RECAL'
include { COMBINE_REPORTS_RECAL } from './Modules/COMBINE_REPORTS_RECAL'
include { MKDUP } from './Modules/MKDUP'
include { CLEANUP_MKDUP_RECAL } from './Modules/CLEANUP_MKDUP_RECAL'
include { CHECK_BAM_MKDUP } from './Modules/CHECK_BAM_MKDUP'
include { COMBINE_REPORTS_MKDUP } from './Modules/COMBINE_REPORTS_MKDUP'
include { MOSDEPTH } from './Modules/MOSDEPTH'
include { CONPAIR } from './Modules/CONPAIR'
include { SAGE } from './Modules/SAGE'
include { STRELKA } from './Modules/STRELKA'
include { MuSE } from './Modules/MuSE'

include { CNVkit_buildcnn } from './Modules/CNVkit_buildcnn'
include { CNVkit } from './Modules/CNVkit'

include { Create_sample_tsv } from './Modules/Delly/Create_sample_tsv'
include { Delly_SVcalling } from './Modules/Delly/Delly_SVcalling'
include { Delly_Prefiltering } from './Modules/Delly/Delly_Prefiltering'
include { Delly_Filtering } from './Modules/Delly/Delly_Filtering'
include { Delly_Filtering_final } from './Modules/Delly/Delly_Filtering_final'

include { MANTA } from './Modules/MANTA'

include { ASCAT } from './Modules/ASCAT'
include { ASCAT_allelecount } from './Modules/ASCAT_exome/ASCAT_allelecount'
include { ASCAT_logrbaf } from './Modules/ASCAT_exome/ASCAT_logrbaf'
include { ASCAT_exome } from './Modules/ASCAT_exome/ASCAT_exome'

include { MUTECT2_CALLING } from './Modules/Mutect2/MUTECT2_CALLING'
include { MUTECT2_CALLING_exome } from './Modules/Mutect2/MUTECT2_CALLING_exome'
include { GETpileUP } from './Modules/Mutect2/GETpileUP'
include { GETpileUP_exome } from './Modules/Mutect2/GETpileUP_exome'
include { GETpileUP_Merge } from './Modules/Mutect2/GETpileUP_Merge'
include { LearnReadOrientationModel } from './Modules/Mutect2/LearnReadOrientationModel'
include { LearnReadOrientationModel_exome } from './Modules/Mutect2/LearnReadOrientationModel_exome'
include { MergeMutectStats } from './Modules/Mutect2/MergeMutectStats'
include { MergeVcfs } from './Modules/Mutect2/MergeVcfs'
include { CalculateContamination } from './Modules/Mutect2/CalculateContamination'
include { FilterMutectCalls } from './Modules/Mutect2/FilterMutectCalls'

include { RENAME_BAM_HEADER } from './Modules/RENAME_BAM_HEADER'
include { GENERATE_INTERVALS } from './Modules/GENERATE_INTERVALS'

include { SUMMARY } from './Modules/SUMMARY.nf'

include { POST } from './Modules/POST.nf'
include { ALLELECOUNTER } from './Modules/ALLELECOUNTER.nf'

include { SAVE_CSV_FASTQC } from './Modules/SAVE_CSV/SAVE_CSV_FASTQC'
include { SAVE_CSV_BWA_MEM } from './Modules/SAVE_CSV/SAVE_CSV_BWA_MEM'
include { SAVE_CSV_MKDUP } from './Modules/SAVE_CSV/SAVE_CSV_MKDUP'
include { SAVE_CSV_RECAL } from './Modules/SAVE_CSV/SAVE_CSV_RECAL'
include { SAVE_CSV_SAGE } from './Modules/SAVE_CSV/SAVE_CSV_SAGE'
include { SAVE_CSV_STRELKA } from './Modules/SAVE_CSV/SAVE_CSV_STRELKA'
include { SAVE_CSV_MuSE } from './Modules/SAVE_CSV/SAVE_CSV_MuSE'
include { SAVE_CSV_Mutect2 } from './Modules/SAVE_CSV/SAVE_CSV_Mutect2'
include { SAVE_CSV_CONPAIR } from './Modules/SAVE_CSV/SAVE_CSV_CONPAIR'
include { SAVE_CSV_MOSDEPTH } from './Modules/SAVE_CSV/SAVE_CSV_MOSDEPTH'

// Define a workflow completion handler
workflow.onComplete {
    def completionStatus = """
==============================================
Pipeline Execution Completed
==============================================
Status        : ${workflow.success ? 'SUCCESS' : 'FAILED'}
Completed at  : ${workflow.complete}
Duration      : ${workflow.duration}
Success       : ${workflow.success}
Exit status   : ${workflow.exitStatus}
Error report  : ${workflow.errorReport ?: 'None'}
==============================================
"""
    writeToLog(params.log_file, completionStatus)
    log.info completionStatus
}

// Helper function to check if a specific tool is selected
def isToolSelected(String tool) {
    if (params.tool) {
        def tools = params.tool.split(',').collect { it.trim() }
        return tools.contains(tool)
    }
    return false
}


workflow {
    chunk = Channel.of(1..20)

    // Check if need to build custom intervals from user provided ref fa file
    if (params.custom_ref) {
    GENERATE_INTERVALS(
        file(params.ref),
        file(params.ref_fai),
        file(params.ref_dict)
    )
    interval_dir_ch = GENERATE_INTERVALS.out.interval_dir
    } else {
        interval_dir_ch = Channel.value(file(params.mutect2_interval_dir))
    }
    
    // Sanity check for ascat input format
    if (params.tool.toString().contains('ascat')) {
        def sampleFile = file(params.sample)
        def hasSexColumn = false
        
        if (sampleFile.exists()) {
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',')
            hasSexColumn = headers.contains("sex")
            
            if (!hasSexColumn) {
                log.error "ERROR: ASCAT tool requires 'sex' column in sample.csv file"
                exit 1
            }
        } else {
            log.error "ERROR: Cannot find input file ${params.sample}"
            exit 1
        }
    }

    // starts from scratch - mapping
    if (params.first_step in ['mapping']) {
        // Sanity check for sample sheet
        def requiredColumns = ['patient', 'sample', 'status', 'fastq_1', 'fastq_2']
        if (params.tool && isToolSelected('ascat')) {
            requiredColumns.add('sex')
        }
        
        // Check if sample sheet exists
        def sampleFile = file(params.sample)
        if (!sampleFile.exists()) {
            log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
            exit 1
        }
        
        // Check for required columns in the sample sheet
        def headerLine = sampleFile.readLines()[0]
        def headers = headerLine.split(',').collect { it.trim() }
        def missingColumns = requiredColumns.findAll { !headers.contains(it) }
        
        if (missingColumns) {
            def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
            log.error errorMsg
            exit 1
        }

        if (params.tool && params.tool.contains('ascat')) {
            Channel.fromPath(params.sample)
            | splitCsv(header:true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'status', 'fastq_1', 'fastq_2', 'sex')
                return meta
            } | set { sample_sheet }
        } else {
            Channel.fromPath(params.sample)
            | splitCsv(header:true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'status', 'fastq_1', 'fastq_2')
                return meta
            } | set { sample_sheet }
        }
        FASTQC(sample_sheet)
        BWA_MEM(sample_sheet).set { BWA_MEM_out }
        CHECK_BAM_BWA(BWA_MEM_out.bam).set{ CHECK_BAM_BWA_out }
        COMBINE_REPORTS_BWA(CHECK_BAM_BWA_out.individual_reports_bwa.collect())
    }

    // starts from markduplicate
    if (params.first_step in ['mapping', 'markdup']) {
        if (params.first_step == "markdup"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet }
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet }
            }
            RENAME_BAM_HEADER(sample_sheet)
            MKDUP(RENAME_BAM_HEADER.out.renamed_bam).set{ MKDUP_out }    
            CLEANUP_BWA(MKDUP.out.cleanup_trigger)
        } else {
            MKDUP(BWA_MEM_out).set{ MKDUP_out }
            CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            CHECK_BAM_MKDUP(MKDUP_out.mkdup_bam).set{ CHECK_BAM_MKDUP_out }
            COMBINE_REPORTS_MKDUP(CHECK_BAM_MKDUP_out.individual_reports_mkdup.collect())
        }

        if (params.genome in ['mm39', 'RN7']){
            MKDUP_out.pair_mutect.filter{it[1].status == 'normal'}.set{normal}
            MKDUP_out.pair_mutect.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            RECALIBRATE_out = [pair_recal: MKDUP_out.pair_mutect]
        }
    }

    // starts from recalibration
    if (params.first_step in ['mapping', 'markdup', 'recalibration']) {
        if (params.first_step == "recalibration" && params.genome in ['GRCh38', 'GRCh37']){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map {row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet_raw }
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet_raw }
            }

            // Rename BAM headers to standardized format
            RENAME_BAM_HEADER(sample_sheet_raw)

            RENAME_BAM_HEADER.out.renamed_bam
            | map { meta, bam, bai ->
                [meta.patient, meta, bam, bai]
            } | set { sample_sheet }

            if (params.type == "exome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(sample_sheet).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_BQSR_out }
                RECALIBRATE_SortBam(RECALIBRATE_BQSR_out.SortBam_input).set{ RECALIBRATE_out }
            } else if (params.type == "genome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal(sample_sheet, chunk, interval_dir_ch).set{ RECALIBRATE_BaseRecal_out }

                RECALIBRATE_BaseRecal_out.MergeReport_input.groupTuple(by:[0,1]).set { BaseRecal_out_pair }
                BaseRecal_out_pair.map{
                    [patient:it[0], status:it[1].status, meta:it[1], bam:it[2][0], table:it[3], bai:it[4][0]]
                }.set{ BaseRecal_out_MAP }

                RECALIBRATE_MergeReport(BaseRecal_out_MAP).set{ RECALIBRATE_MergeReport_out }

                RECALIBRATE_BQSR(RECALIBRATE_MergeReport_out.BQSR_input, chunk, interval_dir_ch).set { RECALIBRATE_BQSR_out }

                RECALIBRATE_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
                BaseRecal_BQSR_out_pair.map{
                    [patient:it[0], meta:it[1], bam:it[2]]
                }.set{BaseRecal_BQSR_out_MAP}

                RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }

                RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }

                CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
            }
        } else {
            if (params.type == "exome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(MKDUP_out.pair_mutect).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_BQSR_out }
                RECALIBRATE_SortBam(RECALIBRATE_BQSR_out.SortBam_input).set{ RECALIBRATE_out }

            } else if (params.type == "genome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal(MKDUP_out.pair_mutect, chunk, interval_dir_ch).set{ RECALIBRATE_BaseRecal_out }

                RECALIBRATE_BaseRecal_out.MergeReport_input.groupTuple(by:[0,1]).set { BaseRecal_out_pair }
                BaseRecal_out_pair.map{
                    [patient:it[0], status:it[1].status, meta:it[1], bam:it[2][0], table:it[3], bai:it[4][0]]
                }.set{ BaseRecal_out_MAP }

                RECALIBRATE_MergeReport(BaseRecal_out_MAP).set{ RECALIBRATE_MergeReport_out }

                RECALIBRATE_BQSR(RECALIBRATE_MergeReport_out.BQSR_input, chunk, interval_dir_ch).set { RECALIBRATE_BQSR_out }

                RECALIBRATE_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
                BaseRecal_BQSR_out_pair.map{
                    [patient:it[0], meta:it[1], bam:it[2]]
                }.set{BaseRecal_BQSR_out_MAP}

                RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }

                RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }

                CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
            }
        }
        CHECK_BAM_RECAL(RECALIBRATE_out.pair_recal).set{ CHECK_BAM_RECAL_out }
        COMBINE_REPORTS_RECAL(CHECK_BAM_RECAL_out.individual_reports_recal.collect())
    }

    // starts from variant calling
    if (params.first_step in ['mapping', 'markdup', 'recalibration', 'variant_calling']) {
        if (params.first_step == "variant_calling"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)  
                | map {row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam), file(row.bai)]
                } | set {sample_sheet_raw}
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet_raw }
            }

            // Rename BAM headers to standardized format
            RENAME_BAM_HEADER(sample_sheet_raw)

            RENAME_BAM_HEADER.out.renamed_bam
            | map { meta, bam, bai ->
                [meta.patient, meta, bam, bai]
            } | set { sample_sheet }

            GETpileUP(sample_sheet, chunk, interval_dir_ch).set { GETpileUP_out }

            sample_sheet.filter{it[1].status == 'normal'}.set{normal}
            sample_sheet.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            // for CN analysis
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], gender:normal[1].sex, normal:normal[2],tumor:tumor[2], sample:tumor[1].sample]
            }.set{ RECALIBRATE_out_MAP_CN }

            // Collect all normal bams for CNVkit and Delly
            normal.map { patient, meta, bam, bai -> [bam, bai] }.flatten().collect().set { normal_bams }

            if (params.genome in ['GRCh38', 'GRCh37']) {
                SAGE(RECALIBRATE_out_MAP)
            }
            
            STRELKA(RECALIBRATE_out_MAP)
            MuSE(RECALIBRATE_out_MAP)
            
            if (params.genome in ['GRCh38', 'GRCh37']) {
                CONPAIR(RECALIBRATE_out_MAP)
                MOSDEPTH(RECALIBRATE_out_MAP)
            }

            MUTECT2_CALLING(RECALIBRATE_out_MAP, chunk, interval_dir_ch).set { MUTECT2_CALLING_out }
            MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
            MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
            MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
            LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
            MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
            MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

            GETpileUP_out.GETpileUP_Merge_input.map{patient, meta, table, chunk -> [patient, meta, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
            GETpileUP_out_pair.map{
            [patient:it[0], meta:it[1], mix:it[2]]
            }.set{ GETpileUP_out_pair_MAP }

            GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
            GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
            }.set{ GETpileUP_out_MAP }

            CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
            def vcfChannel = MergeVcfs_out.MUTECT2_vcf
                .map{ map, vcf ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                }
            def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                .map{ map, table ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, table]]
                }
            def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                .map{ map, segments ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, segments]]
                }
            def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                .map{ map, orientation ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, orientation]]
                }
            def statsChannel = MergeMutectStats_out.MUTECT2_stats
                .map{ map, stats ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                }

            // Join by the composite key (patient_sample)
            vcfChannel
                .join(contaminationChannel)
                .join(segmentsChannel)
                .join(orientationChannel)
                .join(statsChannel)
                .map{ patient_sample, vcfData, contData, segData, oriData, statsData ->
                    [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                }
                .set{ filterInput }

            FilterMutectCalls(filterInput).set{ FILTER_OUT }
            // SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)

            // TOOL-SPECIFIC BLOCKS - MOVED INSIDE variant_calling STEP
            // ASCAT
            if (params.tool && params.tool.contains('ascat')) {
                chromosomes = Channel.of( *(1..22).collect { it.toString() } + ['X'] )
                if (params.type == "exome") {
                    ASCAT_allelecount(RECALIBRATE_out_MAP_CN, chromosomes)

                    ASCAT_allelecount.out.allelecount
                        .map { map, chr, file ->
                            def patient = map.patient
                            def sample = map.sample
                            def gender = map.gender
                            [patient, sample, gender, file]
                        }
                        .groupTuple(by: [0, 1, 2])  // Group by patient, sample, and gender
                        // Restructure from pairs to separate normal and tumor lists
                        .map { patient, sample, gender, files ->
                            // Extract normal and tumor files from the nested structure
                            def normalFiles = []
                            def tumorFiles = []

                            // Flatten any nested structure if present
                            def flatFiles = files.flatten()

                            // Separate normal and tumor files
                            flatFiles.each { file ->
                                if (file.toString().contains("_normal_")) {
                                    normalFiles << file
                                } else if (file.toString().contains("_tumor_")) {
                                    tumorFiles << file
                                }
                            }

                            [patient, sample, gender, normalFiles, tumorFiles]
                        }
                        .set { ascat_logrbaf_input }

                    ASCAT_logrbaf(ascat_logrbaf_input)
                    ASCAT_exome(ASCAT_logrbaf.out.ascat_input)

                } else if (params.type == "genome") {
                    ASCAT(RECALIBRATE_out_MAP_CN)
                }
            }

            // CNVkit
            if (params.tool && params.tool.contains('cnvkit')) {
                CNVkit_buildcnn(normal_bams).set { CNVkit_buildcnn_out }
                CNVkit(RECALIBRATE_out_MAP_CN, CNVkit_buildcnn_out.CNVkit_ref_cnn).set { CNVkit_out }
            }

            // Delly
            if (params.tool && params.tool.contains('delly')) {
                if (params.type == 'exome') {
                    log.error "ERROR: Delly is not recommended for exome data. Please use a different SV caller (manta) for exome data"
                    exit 1  
                } else {
                    Create_sample_tsv(Channel.fromPath(params.sample, checkIfExists: true))
                    sample_tsv_collected = Create_sample_tsv.out.samples_tsv.collect()

                    Delly_SVcalling(RECALIBRATE_out_MAP_CN).set { SVcalling_out }
                    SVcalling_out.SVcalling_bcf
                    .combine(sample_tsv_collected)
                    .set { prefiltering_input }

                    Delly_Prefiltering(prefiltering_input).set { Prefiltering_out }
                    Delly_Filtering(Prefiltering_out.Delly_pre_bcf, normal_bams).set { Filtering_out }

                    Filtering_out.Delly_geno_bcf
                    .combine(sample_tsv_collected)
                    .set { filtering_final_input }

                    Delly_Filtering_final(filtering_final_input)
                }
            }

            // Manta
            if (params.tool && params.tool.contains('manta')) {
                MANTA(RECALIBRATE_out_MAP_CN)
            }

        } else {

            RECALIBRATE_out.pair_recal.filter{it[1].status == 'normal'}.set{normal}
            RECALIBRATE_out.pair_recal.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            // for CN analysis
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], gender:normal[1].sex, normal:normal[2],tumor:tumor[2], sample:tumor[1].sample]
            }.set{ RECALIBRATE_out_MAP_CN }

            // Collect all normal bams for CNVkit and Delly
            normal.map { patient, meta, bam, bai -> [bam, bai] }.flatten().collect().set { normal_bams }

            if (params.genome in ['GRCh38', 'GRCh37']) {
                SAGE(RECALIBRATE_out_MAP)
            }

            STRELKA(RECALIBRATE_out_MAP)
            MuSE(RECALIBRATE_out_MAP)
            
            if (params.genome in ['GRCh38', 'GRCh37']) {
                CONPAIR(RECALIBRATE_out_MAP)
                MOSDEPTH(RECALIBRATE_out_MAP)
            }

            GETpileUP(RECALIBRATE_out.pair_recal, chunk, interval_dir_ch).set { GETpileUP_out }
            MUTECT2_CALLING(RECALIBRATE_out_MAP, chunk, interval_dir_ch).set { MUTECT2_CALLING_out }
            MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
            MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
            MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
            LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
            MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
            MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

            GETpileUP_out.GETpileUP_Merge_input.map{patient, meta, table, chunk -> [patient, meta, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
            GETpileUP_out_pair.map{
            [patient:it[0], meta:it[1], mix:it[2]]
            }.set{ GETpileUP_out_pair_MAP }

            GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
            GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
            }.set{ GETpileUP_out_MAP }

            CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
            def vcfChannel = MergeVcfs_out.MUTECT2_vcf
                .map{ map, vcf ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                }
            def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                .map{ map, table ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, table]]
                }
            def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                .map{ map, segments ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, segments]]
                }
            def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                .map{ map, orientation ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, orientation]]
                }
            def statsChannel = MergeMutectStats_out.MUTECT2_stats
                .map{ map, stats ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                }

            // Join by the composite key (patient_sample)
            vcfChannel
                .join(contaminationChannel)
                .join(segmentsChannel)
                .join(orientationChannel)
                .join(statsChannel)
                .map{ patient_sample, vcfData, contData, segData, oriData, statsData ->
                    [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                }
                .set{ filterInput }

            FilterMutectCalls(filterInput).set{ FILTER_OUT }
            // SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)


            // TOOL-SPECIFIC BLOCKS - ALSO FOR NON variant_calling START STEPS
	    // ASCAT
	    if (params.tool && params.tool.contains('ascat')) {
	        chromosomes = Channel.of( *(1..22).collect { it.toString() } + ['X'] )
	        if (params.type == "exome") {
	            ASCAT_allelecount(RECALIBRATE_out_MAP_CN, chromosomes)
	
	            ASCAT_allelecount.out.allelecount
	                .map { map, chr, file ->
	                    def patient = map.patient
	                    def sample = map.sample
	                    def gender = map.gender
	                    [patient, sample, gender, file]
	                }
	                .groupTuple(by: [0, 1, 2])  // Group by patient, sample, and gender
	                // Restructure from pairs to separate normal and tumor lists
	                .map { patient, sample, gender, files ->
	                    // Extract normal and tumor files from the nested structure
	                    def normalFiles = []
	                    def tumorFiles = []
	
	                    // Flatten any nested structure if present
	                    def flatFiles = files.flatten()
	
	                    // Separate normal and tumor files
	                    flatFiles.each { file ->
	                        if (file.toString().contains("_normal_")) {
	                            normalFiles << file
	                        } else if (file.toString().contains("_tumor_")) {
	                            tumorFiles << file
	                        }
	                    }
	
	                    [patient, sample, gender, normalFiles, tumorFiles]
	                }
	                .set { ascat_logrbaf_input }
	
	            ASCAT_logrbaf(ascat_logrbaf_input)
	            ASCAT_exome(ASCAT_logrbaf.out.ascat_input)
	
	        } else if (params.type == "genome") {
	            ASCAT(RECALIBRATE_out_MAP_CN)
	        }
	    }
	
	    // CNVkit
	    if (params.tool && params.tool.contains('cnvkit')) {
	        CNVkit_buildcnn(normal_bams).set { CNVkit_buildcnn_out }
	        CNVkit(RECALIBRATE_out_MAP_CN, CNVkit_buildcnn_out.CNVkit_ref_cnn).set { CNVkit_out }
	    }
	
	    // Delly
	    if (params.tool && params.tool.contains('delly')) {
	        if (params.type == 'exome' && params.tool.toString().contains('delly')) {
	            log.error "ERROR: Delly is not recommended for exome data. Please use a different SV caller (manta) for exome data"
	            exit 1  
	        } else {
	            Create_sample_tsv(Channel.fromPath(params.sample, checkIfExists: true))
	            sample_tsv_collected = Create_sample_tsv.out.samples_tsv.collect()
	
	            Delly_SVcalling(RECALIBRATE_out_MAP_CN).set { SVcalling_out }
	            SVcalling_out.SVcalling_bcf
	            .combine(sample_tsv_collected)
	            .set { prefiltering_input }
	
	            Delly_Prefiltering(prefiltering_input).set { Prefiltering_out }
	            Delly_Filtering(Prefiltering_out.Delly_pre_bcf, normal_bams).set { Filtering_out }
	
	            Filtering_out.Delly_geno_bcf
	            .combine(sample_tsv_collected)
	            .set { filtering_final_input }
	
	            Delly_Filtering_final(filtering_final_input)
	        }
	    }
	
	    // Manta
	    if (params.tool && params.tool.contains('manta')) {
	        MANTA(RECALIBRATE_out_MAP_CN)
	    }
        }
    }

// POST - filtering and consensus calling
    if (params.first_step in ['mapping', 'markdup', 'recalibration', 'variant_calling'] && 
        params.genome in ['GRCh38', 'GRCh37']) {
        
        // Use the outputs that are guaranteed to exist from variant calling steps
        mutect2Channel = FilterMutectCalls.out.Mutect2_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        museChannel = MuSE.out.MuSE_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        strelkaChannel = STRELKA.out.STRELKA_out
            .map { map, snv, indel ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, snv, indel]]
            }
        
        sageChannel = SAGE.out.SAGE_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        // Recreate recalChannel from the available data
        // Choose the source based on first_step
        sourceChannel = (params.first_step == "variant_calling") ? 
            sample_sheet : 
            RECALIBRATE_out.pair_recal

        // Create recalChannel from the chosen source
        normal_post = sourceChannel.filter{it[1].status == 'normal'}
        tumor_post = sourceChannel.filter{it[1].status == 'tumor'}

        recalChannel = normal_post.cross(tumor_post){it[0]}.map{
            normal, tumor ->
            ["${tumor[1].patient}_${tumor[1].sample}", [tumor[1], tumor[2], tumor[3]]]
        }
        
        postInput = mutect2Channel
            .join(museChannel)
            .join(strelkaChannel)
            .join(sageChannel)
            .join(recalChannel)
            .map { patient_sample, mutect2_vcf, muse_vcf, strelka_vcf, sage_vcf, recal_bam ->
                [patient_sample, mutect2_vcf[1], muse_vcf[1], strelka_vcf[1], strelka_vcf[2], sage_vcf[1], recal_bam[1], recal_bam[2]]
            }
        
        POST(postInput)

        // AlleleCounter for VAF calculation
        def vaf_source = (params.first_step == "variant_calling") ? 
            RENAME_BAM_HEADER.out.renamed_bam.map { meta, bam, bai -> [meta.patient, meta, bam, bai] } :
            RECALIBRATE_out.pair_recal

        vaf_source.filter{it[1].status == 'normal'}.set{normal_vaf}
        vaf_source.filter{it[1].status == 'tumor'}.set{tumor_vaf}

        normal_vaf.cross(tumor_vaf){it[0]}.map{
            normal, tumor ->
            [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1], tumor_bai:tumor[3], normal_bai:normal[3]]
        }.set{ RECALIBRATE_out_MAP_VAF }

        RECALIBRATE_out_MAP_VAF
            .map { map ->
                def patient_sample = "${map.patient}_${map.tumor_meta.sample}"
                [patient_sample, map]
            }
            .set { recal_map_keyed }

        POST.out.final_vcf
            .join(recal_map_keyed)
            .map { patient_sample, vcf_file, recal_map -> 
                [
                    patient_sample: patient_sample,
                    final_vcf: vcf_file,
                    tumor: recal_map.tumor,
                    normal: recal_map.normal,
                    tumor_bai: recal_map.tumor_bai,
                    normal_bai: recal_map.normal_bai
                ]
            }
            .set { allelecounter_input }

        ALLELECOUNTER(allelecounter_input)
    }
}
