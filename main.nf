// Author: Ting Yang, George Wu
// Lab: The Alexandrov Lab @ UCSD
// Date: 2026.9
// Version: 4.0

nextflow.enable.dsl=2

params.database_path = null
params.outdir = "${launchDir}/RESULTS"   
params.log_file = "${launchDir}/pipeline.log"
params.sample = "sample.csv"
params.first_step = null
params.genome = ""
params.tool = ""
params.cleanup = false 
params.mkdup_optical_pixel_distance = 100   // Picard MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE; 2500 for patterned flow cells. Duplicate metrics only.

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

if (!(params.type in ['genome', 'exome'])) {
    error "ERROR: --type must be 'genome' or 'exome' (got: ${params.type == null ? 'not set' : "'${params.type}'"}).\n" +
          "Usage: nextflow run main.nf --type [genome|exome] --genome [GRCh38|GRCh37|mm39|RN7] ..."
}

// --first_step: allowed values, and the genome x first_step combinations that are implemented.
def FIRST_STEPS = ['mapping', 'markdup', 'recalibration', 'variant_calling']
if (!(params.first_step in FIRST_STEPS)) {
    error "ERROR: --first_step must be one of ${FIRST_STEPS.join(', ')} (got: ${params.first_step ?: 'not set'})."
}
// mm39 and RN7 run no base quality score recalibration (their callers use the duplicate-marked BAMs),
// so 'recalibration' is not an entry point for them.
if (params.genome in ['mm39', 'RN7'] && params.first_step == 'recalibration') {
    error "ERROR: --first_step recalibration is not available for genome '${params.genome}': this genome runs no BQSR. " +
          "Start from 'markdup' (raw BAMs) or 'variant_calling' (duplicate-marked BAMs)."
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
            sage_ref_dir: "${params.database_path}/SAGE/hmf_pipeline_resources.38_v2.3.0--2",
            sage_hotspots: "${params.database_path}/SAGE/hmf_pipeline_resources.38_v2.3.0--2/dna/variants/KnownHotspots.somatic.38.vcf.gz",
            sage_panel_bed: "${params.database_path}/SAGE/hmf_pipeline_resources.38_v2.3.0--2/dna/variants/ActionableCodingPanel.38.bed.gz",
            sage_high_confidence_bed: "${params.database_path}/SAGE/hmf_pipeline_resources.38_v2.3.0--2/dna/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz",
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
            sage_ref_dir: "${params.database_path}/SAGE/hmf_pipeline_resources.37_v2.3.0--2",
            sage_hotspots: "${params.database_path}/SAGE/hmf_pipeline_resources.37_v2.3.0--2/dna/variants/KnownHotspots.somatic.37.vcf.gz",
            sage_panel_bed: "${params.database_path}/SAGE/hmf_pipeline_resources.37_v2.3.0--2/dna/variants/ActionableCodingPanel.37.bed.gz",
            sage_high_confidence_bed: "${params.database_path}/SAGE/hmf_pipeline_resources.37_v2.3.0--2/dna/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz",
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
        muse_vcf: "${params.database_path}/Databases/mm39/af_only_mgp_mm39_unique.vcf.gz",
        mutect2_pon: "PoN.mm39.vcf.gz",
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

// Genomes with no validated end-to-end run. RN7 has no Mutect2 panel of normals and has never been
// run through the pipeline: refuse unless the user opts in explicitly.
params.allow_unvalidated_genome = false
if (params.genome == 'RN7' && !params.allow_unvalidated_genome) {
    error "ERROR: genome 'RN7' is not validated: no Mutect2 panel of normals exists and no RN7 sample has been run " +
          "through the pipeline. Re-run with --allow_unvalidated_genome true to proceed; Mutect2 will run without -pon."
}

// Check if user inputs custom ref or not
params.custom_ref = (params.ref && params.ref != params.genomes[params.genome].ref)

if (params.custom_ref) {
    log.warn """
    WARNING: You are using a custom reference genome.
    Ensure your reference uses the same contig naming convention 
    (e.g., chr1 vs 1) as the pipeline's database files for genome '${params.genome}'.
    Mismatched contig names will cause tools like Mutect2, MuSE, BQSR, 
    SAGE, Conpair, and Delly to fail.
    """
}

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
params.mutect2_pon = params.mutect2_pon ?: (params.genomes[params.genome].mutect2_pon ? "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_pon}" : null)
if (params.mutect2_pon && !file(params.mutect2_pon).exists())
    error "ERROR: Mutect2 panel of normals not found: ${params.mutect2_pon}"
if (!params.mutect2_pon)
    log.warn "Genome ${params.genome} has no Mutect2 panel of normals: Mutect2 runs without -pon (lower specificity; not validated for RN7)"
params.mutect2_interval_dir = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_interval_dir}"
params.mutect2_germline = params.mutect2_germline ?: "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_germline}"

params.vaf_script = params.vaf_script ?: "${projectDir}/VAF.py"

// Custom reference: every packaged resource this run will use must share the reference's contig names and lengths (the criterion GATK applies to sequence dictionaries). Any mismatch stops the run before a task is submitted, whatever the consuming tool would do on its own: GATK aborts, but SAGE, mosdepth, MuSE, Strelka, Manta, Delly and ASCAT silently ignore unmatched regions or mis-annotate, which alters the calls.
if (params.custom_ref) {
    def human = params.genome in ['GRCh38', 'GRCh37']
    def tools = params.tool ? params.tool.toString().split(',')*.trim() : []
    def resources = [
        [label: '--ref_dict',       path: params.ref_dict,         kind: 'dict'],
        [label: '--bed',            path: params.bed,              kind: 'bed'],
        [label: 'Mutect2 germline', path: params.mutect2_germline, kind: 'vcf'],
        [label: 'Mutect2 PON',      path: params.mutect2_pon,      kind: 'vcf'],
        [label: 'MuSE dbSNP',       path: params.muse_vcf,         kind: 'vcf'],
    ]
    if (params.type == 'exome')
        resources << [label: 'mosdepth exome BED', path: params.mosdepth_bed, kind: 'bed']
    if (human && params.first_step != 'variant_calling') {
        resources << [label: 'BQSR known sites 1', path: params.recal_knownsite1, kind: 'vcf']
        resources << [label: 'BQSR known sites 2', path: params.recal_knownsite2, kind: 'vcf']
        if (params.type == 'exome')
            resources << [label: 'BQSR exome intervals', path: params.recal_interval_wes, kind: 'interval_list']
    }
    if (human) {
        resources << [label: 'SAGE hotspots',            path: params.SAGE_hotspots,            kind: 'vcf']
        resources << [label: 'SAGE panel BED',           path: params.SAGE_panel_bed,           kind: 'bed']
        resources << [label: 'SAGE high-confidence BED', path: params.SAGE_high_confidence_bed, kind: 'bed']
        resources << [label: 'Conpair markers',          path: params.conpair_marker,           kind: 'bed']
    }
    if ('delly' in tools)
        resources << [label: 'Delly exclude list', path: params.delly_excl, kind: 'bed', allowUnknown: true]
    if ('ascat' in tools && human) {
        // alleleCounter loci (chromosome, position per line), one file per chromosome. The allele, GC-content
        // and replication-timing files of the same bundle share these coordinates.
        def loci_prefix = params.type == 'exome' ? params.ascat_loci_prefix_wes : params.ascat_loci_prefix
        ((1..22).collect { it.toString() } + ['X']).each { c ->
            resources << [label: "ASCAT loci chr${c}", path: "${loci_prefix}${c}.txt", kind: 'bed', maxLines: 5000]
        }
    }

    def check = ReferenceCheck.run(params.ref_fai, resources)
    def problems = check.errors + check.warnings
    if (problems) {
        error "ERROR: custom reference ${params.ref} is not compatible with the packaged resources for genome '${params.genome}':\n  - " +
              problems.join('\n  - ') +
              "\nProvide resources built for this reference through the matching --<resource> options, " +
              "leave the tool out of --tool, or use the packaged reference."
    }
    log.info "Custom reference check: ${resources.count { it.path }} packaged resources checked against ${params.ref}: all contig names and lengths match"
}

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
Output Dir     : ${params.outdir}
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

params.bam_dir="${params.outdir}/BAM"
params.report_dir="${params.outdir}/REPORT"
params.manta_bed="${params.database_path}/GRCh38_ref/manta.bed.gz"
params.mkdup_temp_dir="${launchDir}/mkdup_tmp"
params.FASTQC_dir="${params.outdir}/FASTQC"
params.mkdup_dir="${params.outdir}/MKDUP"
params.recal_dir="${params.outdir}/RECALIBRATE"
params.SAGE_jar="${params.database_path}/SAGE/sage_v4.2.jar"
params.SAGE_dir="${params.outdir}/SAGE"
params.strelka_dir="${params.outdir}/STRELKA"
params.muse_dir="${params.outdir}/MuSE"
params.MUTECT2_dir="${params.outdir}/Mutect2"
params.mosdepth_dir="${params.outdir}/mosdepth"
params.conpair_dir="${params.outdir}/Conpair"
params.cnvkit_dir="${params.outdir}/CNVkit"
params.Delly_dir="${params.outdir}/Delly"
params.ascat_dir="${params.outdir}/ASCAT"
params.conpair="${params.database_path}/Conpair-0.2"
params.manta_dir="${params.outdir}/MANTA"
params.post_dir="${params.outdir}/POST"
params.intervals_dir="${params.outdir}/custom_intervals"
params.VAF_dir="${params.outdir}/POST/VAF"

params.bwamem2_env = "$projectDir/yml/bwamem2.yml"
params.mkdup_env = "$projectDir/yml/picard.yml"
params.conpair_env = "$projectDir/yml/conpair.yml"
params.samtools_env = "$projectDir/yml/samtools.yml"
params.strelka_env = "$projectDir/yml/strelka_env.yml"
params.mosdepth_env = "$projectDir/yml/mosdepth_env.yml"
params.fastqc_env = "$projectDir/yml/fastqc_env.yml"
params.cnvkit_env = "$projectDir/yml/cnvkit.yml"
params.delly_env = "$projectDir/yml/delly.yml"
params.ascat_env = "$projectDir/yml/ascat.yml"
params.alleleCounter_env = "$projectDir/yml/alleleCounter.yml"
params.manta_env = "$projectDir/yml/manta.yml"
params.SAGE_java_env = "$projectDir/yml/sage_java.yml"
params.muse_env = "$projectDir/yml/muse.yml"
params.picard_merge_env = "$projectDir/yml/picard.yml"
params.gatk_env = "${projectDir}/yml/gatk_runtime.yml"
params.POST_env = "${projectDir}/yml/POST.yml"
params.gatk = "${params.database_path}/gatk-4.6.2.0/gatk"

// Write configuration to log file
writeToLog(params.log_file, configInfo)

// Also display in terminal
log.info configInfo

if (params.cleanup) {
    def cleanupWarning = """
==============================================
WARNING: --cleanup is ENABLED
==============================================
Intermediate files will be permanently deleted during the run:
  - BWA_MEM work directories and raw BAMs
  - MKDUP BAMs and work directories
  - Recalibration chunk/merged BAMs and work directories

You will NOT be able to use -resume to reuse these steps.
Only final outputs in RESULTS/ will be kept.
==============================================
"""
    log.warn cleanupWarning
    writeToLog(params.log_file, cleanupWarning)
}

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
include { GETpileUP } from './Modules/Mutect2/GETpileUP'
include { GETpileUP_Merge } from './Modules/Mutect2/GETpileUP_Merge'
include { LearnReadOrientationModel } from './Modules/Mutect2/LearnReadOrientationModel'
include { MergeMutectStats } from './Modules/Mutect2/MergeMutectStats'
include { MergeVcfs } from './Modules/Mutect2/MergeVcfs'
include { CalculateContamination } from './Modules/Mutect2/CalculateContamination'
include { FilterMutectCalls } from './Modules/Mutect2/FilterMutectCalls'

include { RENAME_BAM_HEADER } from './Modules/RENAME_BAM_HEADER'
include { GENERATE_INTERVALS } from './Modules/GENERATE_INTERVALS'
include { POST_SCENARIOS } from './Modules/POST_scenarios.nf'
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
    def ignored = workflow.stats.ignoredCount
    def status  = !workflow.success ? 'FAILED'
                : ignored > 0       ? "COMPLETED WITH ${ignored} IGNORED TASK FAILURE(S) - outputs missing, see .nextflow.log"
                :                     'SUCCESS'
    def completionStatus = """
==============================================
Pipeline Execution Completed
==============================================
Status        : ${status}
Tasks         : ${workflow.stats.succeedCount} succeeded, ${workflow.stats.failedCount} failed, ${ignored} ignored
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

def validatePairing(sheet) {
    def sampleFile = file(sheet)
    if (!sampleFile.exists()) error "ERROR: Sample sheet file '${sheet}' does not exist"
    def header = sampleFile.readLines()[0].split(',')*.trim()
    if (!['patient', 'sample', 'status'].every { it in header }) return      // reported by the per-step column check
    def rows = sampleFile.splitCsv(header: true)
    def problems = []
    // patient and sample become file names, read-group fields and shell arguments: letters, digits, '.', '_' and '-' only, starting with a letter or digit (a leading '-' could be read as an option). Values are checked untrimmed so that stray spaces are rejected rather than carried into names.
    def ID_OK = ~/^[A-Za-z0-9][A-Za-z0-9._-]*$/
    rows.eachWithIndex { r, i ->
        ['patient', 'sample'].each { col ->
            def v = r[col]
            if (!v || !(v ==~ ID_OK))
                problems << "row ${i + 2}: ${col} '${v}' is not a valid identifier (letters, digits, '.', '_', '-' only, starting with a letter or digit)"
        }
    }
    rows.groupBy { [it.patient?.trim(), it.sample?.trim(), it.status?.trim()] }.each { k, v ->
        if (v.size() > 1) problems << "patient ${k[0]}, sample ${k[1]}, status ${k[2]}: listed ${v.size()} times. " +
            "One FASTQ pair (or one BAM) per sample: concatenate the lanes of a multi-lane library into one pair first (README, 'One FASTQ pair per sample')"
    }
    rows.groupBy { it.patient?.trim() }.each { patient, prows ->
        def status  = prows.collect { it.status }
        def unknown = status.findAll { !(it in ['tumor', 'normal']) }.unique()
        if (unknown)                  problems << "patient ${patient}: status must be 'tumor' or 'normal' (got: ${unknown.join(', ')})"
        def nNormal = status.count('normal'), nTumor = status.count('tumor')
        if (nNormal != 1)             problems << "patient ${patient}: ${nNormal} normal sample(s), exactly 1 is required"
        if (nTumor == 0)              problems << "patient ${patient}: no tumor sample"
    }
    if (problems) {
        error "ERROR: sample sheet ${sheet} is not valid:\n  - " + problems.join('\n  - ') +
              "\nEach patient needs exactly one row with status 'normal' and one or more with status 'tumor'."
    }
    log.info "Sample sheet ${sheet}: ${rows.count { it.status?.trim() == 'tumor' }} tumor/normal pair(s) across ${rows*.patient.unique().size()} patient(s)"
}


workflow {
    validatePairing(params.sample)
    chunk = Channel.of(1..20)

    // BAMs removed by samtools quickcheck (CHECK_BAM_BWA / CHECK_BAM_MKDUP), as
    // tuple(patient, meta, stage). Used by the pair filter after CHECK_BAM_RECAL.
    qc_failed = Channel.empty()

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

        Channel.fromPath(params.sample)
        | splitCsv(header:true)
        | map { row ->
            def keys = ['patient', 'sample', 'status'] + (isToolSelected('ascat') ? ['sex'] : [])
            def meta = row.subMap(keys)
            [meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
        } | set { sample_sheet }
        FASTQC(sample_sheet)
        BWA_MEM(sample_sheet).set { BWA_MEM_out }
        CHECK_BAM_BWA(BWA_MEM_out.bam).set{ CHECK_BAM_BWA_out }
        COMBINE_REPORTS_BWA(CHECK_BAM_BWA_out.individual_reports_bwa.collect())

        // Number of BAMs per patient in the sample sheet (used by the pair filter)
        sample_sheet
            .map { meta, fq1, fq2 -> meta.patient }
            .toList()
            .map { patients -> patients.countBy { it } }
            .set { bams_per_patient }

        // Only BAMs that pass quickcheck go on to MarkDuplicates
        CHECK_BAM_BWA_out.checked
            .branch { meta, bam, qc ->
                pass: qc == 'PASS'
                fail: true
            }
            .set { bwa_qc }
        bwa_qc.pass
            .map { meta, bam, qc -> tuple(meta, bam) }
            .set { BWA_MEM_good }
        qc_failed = qc_failed.mix(
            bwa_qc.fail.map { meta, bam, qc ->
                log.warn "quickcheck FAIL after BWA: removing ${meta.patient}_${meta.sample}_${meta.status} (${bam.name})"
                tuple(meta.patient, meta, 'BWA')
            }
        )
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
            sample_sheet
                .map { meta, bam, bai -> meta.patient }
                .toList()
                .map { patients -> patients.countBy { it } }
                .set { bams_per_patient }

            RENAME_BAM_HEADER(sample_sheet)
            MKDUP(RENAME_BAM_HEADER.out.renamed_bam).set{ MKDUP_out }    
            if (params.cleanup) {
                CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            }
        } else {
            MKDUP(BWA_MEM_good).set{ MKDUP_out }
            if (params.cleanup) {
                CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            }            
        }

        CHECK_BAM_MKDUP(MKDUP_out.mkdup_bam).set{ CHECK_BAM_MKDUP_out }
        COMBINE_REPORTS_MKDUP(CHECK_BAM_MKDUP_out.individual_reports_mkdup.collect())

        // Only BAMs that pass quickcheck go on to recalibration.
        // MKDUP_good has the same shape as MKDUP's pair_mutect output: tuple(patient, meta, bam, bai)
        CHECK_BAM_MKDUP_out.checked
            .branch { meta, bam, bai, qc ->
                pass: qc == 'PASS'
                fail: true
            }
            .set { mkdup_qc }
        mkdup_qc.pass
            .map { meta, bam, bai, qc -> tuple(meta.patient, meta, bam, bai) }
            .set { MKDUP_good }
        qc_failed = qc_failed.mix(
            mkdup_qc.fail.map { meta, bam, bai, qc ->
                log.warn "quickcheck FAIL after MarkDuplicates: removing ${meta.patient}_${meta.sample}_${meta.status} (${bam.name})"
                tuple(meta.patient, meta, 'MKDUP')
            }
        )

        if (params.genome in ['mm39', 'RN7']){
            MKDUP_good.filter{it[1].status == 'normal'}.set{normal}
            MKDUP_good.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            RECALIBRATE_out = [pair_recal: MKDUP_good]
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

            sample_sheet_raw
                .map { meta, bam, bai -> meta.patient }
                .toList()
                .map { patients -> patients.countBy { it } }
                .set { bams_per_patient }

            // Rename BAM headers to standardized format
            RENAME_BAM_HEADER(sample_sheet_raw)

            RENAME_BAM_HEADER.out.renamed_bam
            | map { meta, bam, bai ->
                [meta.patient, meta, bam, bai]
            } | set { sample_sheet }

            if (params.type == "exome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(sample_sheet).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_out }
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

                RECALIBRATE_MergeBam_out.MergeBam_input
                    .map { patient, meta, merged -> [[patient, meta.sample, meta.status], patient, meta, merged] }
                    .join( RECALIBRATE_MergeReport_out.BQSR_input.map { patient, meta, bam, table, bai -> [[patient, meta.sample, meta.status], bam, bai] }, failOnMismatch: true )
                    .map { key, patient, meta, merged, mkdup_bam, mkdup_bai -> [patient, meta, merged, mkdup_bam, mkdup_bai] }
                    .set { SortBam_input }
                RECALIBRATE_SortBam(SortBam_input).set{ RECALIBRATE_out }

                if (params.cleanup) {
                    CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
                }
            }
        } else {
            if (params.type == "exome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(MKDUP_good).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_out }
            } else if (params.type == "genome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal(MKDUP_good, chunk, interval_dir_ch).set{ RECALIBRATE_BaseRecal_out }

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

                RECALIBRATE_MergeBam_out.MergeBam_input
                    .map { patient, meta, merged -> [[patient, meta.sample, meta.status], patient, meta, merged] }
                    .join( RECALIBRATE_MergeReport_out.BQSR_input.map { patient, meta, bam, table, bai -> [[patient, meta.sample, meta.status], bam, bai] }, failOnMismatch: true )
                    .map { key, patient, meta, merged, mkdup_bam, mkdup_bai -> [patient, meta, merged, mkdup_bam, mkdup_bai] }
                    .set { SortBam_input }
                RECALIBRATE_SortBam(SortBam_input).set{ RECALIBRATE_out }

                if (params.cleanup) {
                    CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
                }
            }
        }
        CHECK_BAM_RECAL(RECALIBRATE_out.pair_recal).set{ CHECK_BAM_RECAL_out }
        COMBINE_REPORTS_RECAL(CHECK_BAM_RECAL_out.individual_reports_recal.collect())

        // Pair filter: keep only BAMs that passed every quickcheck, and only for
        // patients that still have at least one normal and one tumor BAM.
        // Every BAM in the sample sheet ends here exactly once, either as a
        // recalibrated BAM or as a quickcheck failure from an earlier step, so each
        // patient is released as soon as all of its own BAMs are accounted for.
        CHECK_BAM_RECAL_out.checked
            .map { patient, meta, bam, bai, qc ->
                if (qc != 'PASS') {
                    log.warn "quickcheck FAIL after recalibration: removing ${meta.patient}_${meta.sample}_${meta.status} (${bam.name})"
                }
                tuple(patient, meta, qc == 'PASS' ? [bam, bai] : null, qc == 'PASS' ? 'PASS' : 'RECAL')
            }
            .mix(qc_failed.map { patient, meta, stage -> tuple(patient, meta, null, stage) })
            .combine(bams_per_patient)
            .map { patient, meta, files, stage, counts -> tuple(groupKey(patient, counts[patient]), meta, files, stage) }
            .groupTuple()
            .flatMap { patient, metas, files, stages ->
                def rows = [metas, files, stages].transpose()
                def good = rows.findAll { m, f, s -> s == 'PASS' }
                def has_normal = good.any { m, f, s -> m.status == 'normal' }
                def has_tumor  = good.any { m, f, s -> m.status == 'tumor' }
                if (!(has_normal && has_tumor)) {
                    def missing = [has_normal ? null : 'normal', has_tumor ? null : 'tumor'].findAll()
                    log.warn "Skipping patient ${patient}: no passing ${missing.join(' or ')} BAM"
                    return []
                }
                good.collect { m, f, s -> tuple(m.patient, m, f[0], f[1]) }
            }
            .set { RECAL_good }
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
                    .map { map, chr, tumor_file, normal_file -> [map.patient, map.sample, map.gender, normal_file, tumor_file] }
                    .groupTuple(by: [0, 1, 2])   // -> [patient, sample, gender, normalFiles, tumorFiles]
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

            RECAL_good.filter{it[1].status == 'normal'}.set{normal}
            RECAL_good.filter{it[1].status == 'tumor'}.set{tumor}
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

            GETpileUP(RECAL_good, chunk, interval_dir_ch).set { GETpileUP_out }
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
                    .map { map, chr, tumor_file, normal_file -> [map.patient, map.sample, map.gender, normal_file, tumor_file] }
                    .groupTuple(by: [0, 1, 2])   // -> [patient, sample, gender, normalFiles, tumorFiles]
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
            RECAL_good

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
		
        POST_SCENARIOS(postInput)
        
        // AlleleCounter for VAF calculation
        def vaf_source = (params.first_step == "variant_calling") ? 
            RENAME_BAM_HEADER.out.renamed_bam.map { meta, bam, bai -> [meta.patient, meta, bam, bai] } :
            RECAL_good

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

        POST_SCENARIOS.out.final_vcf_s4
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
            .set { allelecounter_input_4 }


        POST_SCENARIOS.out.final_vcf_s5
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
            .set { allelecounter_input_5 }

        POST_SCENARIOS.out.final_vcf_s6
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
            .set { allelecounter_input_6 }


        ch_out1 = allelecounter_input_4.map { file -> tuple(file, 'scenario_4') }
        ch_out2 = allelecounter_input_5.map { file -> tuple(file, 'scenario_5') }
        ch_out3 = allelecounter_input_6.map { file -> tuple(file, 'scenario_6') }

        ch_all_outputs = ch_out1.mix(ch_out2, ch_out3)

        ALLELECOUNTER(ch_all_outputs)
    }
}
