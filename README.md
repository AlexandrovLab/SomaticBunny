# SomaticBunny

SomaticBunny is an ensemble pipeline for somatic variant calling that integrates multiple tools to provide comprehensive variant detection in tumor-normal paired samples.

<p align="center">
  <img src="https://github.com/AlexandrovLab/SomaticBunny/blob/main/workflow_logo/SomaticBunny.png" alt="SomaticBunny Pipeline Workflow"/>
</p>

## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Basic Usage](#basic-usage)
  - [Custom Reference Genome](#custom-reference-genome)
  - [Advanced Options](#advanced-options)
- [Output](#output)
- [Tool Versions](#tool-versions)

## Overview

The SomaticBunny pipeline integrates multiple variant calling tools to improve detection accuracy for:
- Single Nucleotide Variants (SNVs)
- Insertions and Deletions (INDELs)
- Structural Variants (SVs)
- Copy Number Variants (CNVs)

This pipeline is designed to work with paired tumor-normal whole-genome or whole-exome sequencing data.

## Prerequisites

- Conda/Mamba for environment management

## Installation

1. Clone this repository

2. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment

3. Make sure you have the following required files in your working directory:
   - `main.nf`
   - `nextflow.config`
   - `conf/base.config`
   - `sample.csv`

4. Prepare your `sample.csv` file according to your analysis requirements:

   <details>
   <summary><b>For starting from fastq (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,fastq_1,fastq_2
   UCSD101,A,normal,/PATH/TO/ERR5285401_1.fastq.gz,/PATH/TO/ERR5285401_2.fastq.gz,XY
   UCSD101,A,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz,XY
   UCSD101,B,tumor,/PATH/TO/ERR5285403_1.fastq.gz,/PATH/TO/ERR5285403_2.fastq.gz,XY
   ```
   </details>
   <details>
   <summary><b>For starting from bams (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,bam,bai,sex
   UCSD101,A,normal,/PATH/TO/UCSD101_A_normal.bam,/PATH/TO/UCSD101_A_normal.bam.bai,XY
   UCSD101,A,tumor,/PATH/TO/UCSD101_A_tumor.bam,/PATH/TO/UCSD101_A_tumor.bam.bai,XY
   UCSD101,B,tumor,/PATH/TO/UCSD101_B_tumor.bam,/PATH/TO/UCSD101_B_tumor.bam.bai,XY
   ```
   </details>

   <details>
   <summary><b>For ASCAT analysis (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,fastq_1,fastq_2,sex
   UCSD101,A,normal,/PATH/TO/ERR5285401_1.fastq.gz,/PATH/TO/ERR5285401_2.fastq.gz,XY
   UCSD101,A,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz,XY
   UCSD101,B,tumor,/PATH/TO/ERR5285403_1.fastq.gz,/PATH/TO/ERR5285403_2.fastq.gz,XY
   ```
   </details>



## Configuration

For memory management, modify the following parameters if needed:

- Set `$params.mkdup_temp_dir` in `main.nf` to your ideal location (Default: `$projectDir/mkdup_tmp`)
- Set `$workDir` in `nextflow.config` to your ideal location (Default: `./work`)

```nextflow
// In main.nf
params.mkdup_temp_dir = "/path/to/ideal/mkdup_tmp/folder/"

// In nextflow.config
workDir = "/path/to/ideal/intermediate/work/folder/"
```

## Usage

### Basic Usage

1. Request an interactive compute node:

```bash
srun -N 1 -n 1 -c 16 --mem 250G -t 80:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash
```

2. Activate your Nextflow environment:

```bash
conda activate env_nf
```

3. Set TSCC temporary directory and Nextflow tmp directory, and memory:

```bash
export TMPDIR=/path/to/restricted/folder/

export NXF_OPTS="-Djava.io.tmpdir=${TMPDIR} -Xms4g -Xmx16g"
```

4. Run the pipeline:

```bash
nextflow run main.nf --type [genome|exome] --genome [GRCh38|GRCh37|mm39|RN7] --first_step [mapping|markdup|recalibration|variant_calling] --tool [ascat,manta,...] --database_path [/path/to/SomaticBunny_database] --ref [/path/to/ref/] --ref_fai [/path/to/ref_fai] --ref_dict [/path/to/ref_dict] --bed [/path/to/bed] --bed_tbi [/path/to/bed_tbi]
```

#### Pipeline Steps and Tools

The SomaticBunny pipeline workflow is divided into multiple steps that can be run individually or in sequence using the `--first_step` parameter.

##### Available first_step

| First_step | Description | Required Input | 
|------|-------------| ---------------|
| `mapping` | This will start from performing alignment of FASTQ files to the reference genome | fastq_1 and fastq_2 |
| `markdup` | This will start from mark duplicate reads using the aligned BAM files | raw bam and bai |
| `recalibration` | This will start from performing base quality score recalibration | mark duplicated bam and bai |
| `variant_calling` | This will start from executing variant calling tools | recalibrated bam and bai |

##### Available reference genomes
| Genomes | 
|---------|
| GRCh38  | 
| GRCh37  | 
| mm39  | 
| RN7  | 

##### Available Tools

Use the `--tool` parameter to specify which variant callers to run. You can select multiple tools by separating them with commas.

| Category | Available Tools | Description |
|----------|----------------|-------------|
| SV Callers | `manta` | Structural variant and indel caller |
| | `delly` | Integrated structural variant detection |
| CNV Callers | `ascat` | Allele-specific copy number analysis (requires sex information) |
| | `cnvkit` | Copy number variation detection from targeted DNA sequencing |

Example of running multiple tools:
```bash
# Run Manta for SVs and ASCAT for CNVs
nextflow run main.nf --genome GRCh38 --type genome --first_step mapping --tool manta,ascat
```

### Custom Reference Genome

When starting the pipeline from BAM files (`--first_step markdup`, `recalibration`, or `variant_calling`), you **must** provide your own reference genome files. This is because BAM files are already aligned to a specific reference, and the pipeline needs the matching reference for downstream analysis.

#### Required parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--ref` | Reference FASTA file | `/path/to/hg38_genome.fa` |
| `--ref_fai` | FASTA index file | `/path/to/hg38_genome.fa.fai` |
| `--ref_dict` | Sequence dictionary | `/path/to/hg38_genome.dict` |
| `--bed` | Callable regions BED file (bgzipped) | `/path/to/callable_regions.bed.gz` |
| `--bed_tbi` | Tabix index for the BED file | `/path/to/callable_regions.bed.gz.tbi` |

All five files must be provided together. The pipeline will exit with an error if any are missing.

#### Example

```bash
nextflow run main.nf \
  --type genome \
  --genome GRCh38 \
  --first_step variant_calling \
  --sample sample.csv \
  --database_path /path/to/SomaticBunny_database \
  --ref /path/to/hg38_genome.fa \
  --ref_fai /path/to/hg38_genome.fa.fai \
  --ref_dict /path/to/hg38_genome.dict \
  --bed /path/to/hg38_genome_callable_regions.bed.gz \
  --bed_tbi /path/to/hg38_genome_callable_regions.bed.gz.tbi
```

When a custom reference is provided, the pipeline will automatically generate the scattered interval lists required by GATK tools (Mutect2, BaseRecalibrator, etc.) from your reference. When starting from `mapping`, these parameters are optional and will fall back to the pipeline's built-in defaults for the selected `--genome`.

> [!IMPORTANT]
> **Your custom reference genome must use the `chr` prefix for contig names** (e.g., `chr1`, `chr2`, ..., `chrX`, `chrY`). The pipeline's built-in database files (dbSNP, gnomAD, Panel of Normals, known indels, etc.) use `chr`-prefixed contig names. A reference without the `chr` prefix (e.g., `1`, `2`, ..., `X`, `Y`) will cause tools such as Mutect2, MuSE, Strelka, SAGE, and BQSR to fail due to contig name mismatches.
>
> Different GRCh38 builds that include additional or fewer contigs (e.g., alt contigs, decoys, HLA) are supported as long as the standard chromosomes use the `chr` prefix.

### Advanced Options

#### Publishing Intermediate Files

By default, intermediate files are not saved to reduce disk usage. Use the `--publish` parameter to save specific file types:

```bash
nextflow run main.nf --type genome --genome GRCh38 --first_step mapping --publish raw_bam,markdup_bam
```

**Available options:**
- `raw_bam`: Publish raw BAM files from BWA-MEM
- `markdup_bam`: Publish mark duplicates BAM files
- `recal_bam`: Publish recalibrated BAM files

#### Resuming Failed Runs

If your pipeline terminates with an error or the interactive node is killed, resume with:

```bash
nextflow run main.nf --type genome -genome GRCh38 --first_step mapping --tool ascat,manta -resume
```

#### Email Notifications

Receive completion notification:

```bash
nextflow run main.nf --type genome -genome GRCh38 --first_step mapping --tool ascat,manta -N your_email@example.com
```

## Output

All results and reports are stored in the **RESULTS** folder with the following structure:

```
RESULTS/
├── ASCAT (optional)
├── BAM (optional)
├── CNVkit (optional)
├── Conpair
├── Delly (optional)
├── FASTQC
├── MANTA (optional)
├── MKDUP
├── mosdepth
├── MuSE
├── Mutect2
├── RECALIBRATE (optional)
├── REPORT
├── SAGE
└── STRELKA
```

## Tool Versions

Primary tools and compatibility-critical runtime dependencies are
version-pinned in the Conda environment files under `yml/`. External
software distributed with the database bundle is documented separately
by SHA256 checksum.

| Category                    | Tool                                         | Version  |
| --------------------------- | -------------------------------------------- | -------- |
| Quality Control             | FastQC                                       | v0.12.1  |
| Read Alignment              | bwa-mem2                                     | v2.2.1   |
| BAM/VCF Processing          | samtools                                     | v1.21    |
|                             | HTSlib                                       | v1.21    |
|                             | Picard AddOrReplaceReadGroups/MarkDuplicates | v3.2.0   |
|                             | mosdepth                                     | v0.3.8   |
| Sample Identity/QC          | Conpair                                      | v0.2     |
|                             | Python (Conpair runtime)                     | v2.7.15  |
|                             | NumPy (Conpair dependency)                   | v1.16.5  |
|                             | SciPy (Conpair dependency)                   | v1.2.1   |
|                             | Matplotlib (Conpair dependency)              | v2.2.5   |
|                             | pandas (Conpair dependency)                  | v0.24.2  |
|                             | pysam (Conpair dependency)                   | v0.15.4  |
|                             | GATK (Conpair dependency)                    | v3.8     |
|                             | OpenJDK (Conpair runtime)                    | v8.0.412 |
|                             | samtools (Conpair dependency)                | v1.9     |
|                             | HTSlib (Conpair dependency)                  | v1.9     |
| SNV/INDEL Calling           | Strelka2                                     | v2.9.10  |
|                             | Mutect2 (GATK4)                              | v4.6.0.0 |
|                             | SAGE                                         | v3.3     |
|                             | MuSE2                                        | v2.1.2   |
| GATK4 Runtime               | OpenJDK                                      | v17.0.18 |
|                             | Python                                       | v3.10.13 |
|                             | HTSJDK (bundled with GATK)                   | v4.1.1   |
|                             | Picard (bundled with GATK)                   | v3.2.0   |
| Mutect2 VCF Processing      | Picard MergeVcfs                             | v2.18.27 |
| SAGE Runtime                | OpenJDK                                      | v17.0.18 |
| Variant Filtering           | DKFZ bias filter                             | v1.2.3a  |
| Structural Variant Calling  | Delly                                        | v1.3.1   |
|                             | BCFtools                                     | v1.17    |
|                             | Manta                                        | v1.6.0   |
| Copy Number Analysis        | CNVkit                                       | v0.9.8   |
|                             | Python (CNVkit runtime)                      | v3.9.20  |
|                             | NumPy (CNVkit dependency)                    | v1.23.5  |
|                             | pandas (CNVkit dependency)                   | v1.5.3   |
|                             | Matplotlib (CNVkit dependency)               | v3.7.3   |
|                             | R (CNVkit runtime)                           | v4.3.3   |
|                             | DNAcopy                                      | v1.76.0  |
|                             | ASCAT                                        | v3.1.1   |
|                             | R (ASCAT runtime)                            | v4.2.2   |
|                             | cancerit-allelecount                         | v4.3.0   |
| Python Variant Processing   | pysam                                        | v0.24.0  |
| Reporting and Summarization | Python                                       | v3.13.0  |
|                             | NumPy                                        | v2.1.2   |
|                             | pandas                                       | v2.2.3   |

## Software Provenance

The following external software artifacts are distributed with the
database bundle rather than installed directly by the Conda environment
files.

| Software      | Version  | Relative artifact                                  | SHA256                                                             |
| ------------- | -------- | -------------------------------------------------- | ------------------------------------------------------------------ |
| Conpair       | v0.2     | `Conpair-0.2/` source tree                         | `d144e413959d05a191375a14925e8dc6c90e4ab4fba514456819b0472dce5602` |
| SAGE          | v3.3     | `SAGE/sage_v3.3.jar`                               | `290df227c91b15ef3c932cf66b2c7cd52c617071233a7957ca7e6e54889e2b79` |
| GATK launcher | v4.6.0.0 | `gatk-4.6.0.0/gatk`                                | `91fa870b65312522f8c77f0f7adfb5a92b8da97d91c05302ce3535516036d1b8` |
| GATK package  | v4.6.0.0 | `gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar`      | `8b7d8b78ad1ac8f916176e717ea5bc649576080208f608795f53f6fa264bf0ee` |

The Conpair checksum was calculated from the sorted SHA256 manifest of
all source files, excluding `.git`, `__pycache__`, and Python bytecode
files. File-level checksums for SAGE and GATK were calculated using
`sha256sum`.

## Database Provenance

Reference genomes and supporting database resources are supplied
separately through the `--database_path` parameter. They are not
downloaded or installed while the pipeline is running.

All paths below are relative to the directory supplied through
`--database_path`. SHA256 checksums for all 628 database files,
including reference indexes, VCF indexes, interval shards, ASCAT
resources, and SAGE reference files, are recorded in
[`checksums/database_files.sha256`](checksums/database_files.sha256).

### Primary Reference and Database Files

| Category | Resource | Relative artifact | SHA256 |
| -------- | -------- | ----------------- | ------ |
| Reference genome | GRCh38 (`GRCh38.d1.vd1`) FASTA | `GRCh38_ref/GRCh38.d1.vd1.fa` | `fbde121e898f00cc064eaafe5acc58940d7fe069ad14d34acd40c5d675e6b69b` |
| Genomic intervals | GRCh38 chromosome BED | `GRCh38_ref/hg38_chr.bed.gz` | `e9a71588bcc3a556296c9b8f200143e54b05a889cded3aaf4fcd3384a1f0c470` |
| Sample identity/QC | Conpair 1000 Genomes Phase 3 markers, GRCh38 liftover, 20130502, MAF 0.4, LD 0.8 | `Databases/GRCh38/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed` | `8dff7382c5eb59b498a39c40f2925c742c8e4270a70e055275607d665dc2d8d2` |
| Sample identity/QC | Conpair 1000 Genomes Phase 3 genotypes, GRCh38 liftover, 20130502, MAF 0.4, LD 0.8 | `Databases/GRCh38/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt` | `1e9d9c21e0ed78a55cd98894679ab80e5de6fdf138610599a10afca6ad423114` |
| Coverage analysis | GRCh38 exome target BED | `Databases/GRCh38/GRCh38_exome.bed` | `e7032d588a9dffe74f3c6f31cb0b5de31f2b4e20b8a78ec1024acece1edb2ae6` |
| Known variants | dbSNP build 138, GRCh38 | `Databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz` | `f406ebcc1b877b17cd3ad56a63b4e570e5614524e495e9614d56ecfb05de19eb` |
| BQSR known sites | Broad hg38 dbSNP build 138 | `Databases/GRCh38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf` | `e9a8c26560b7c1afe82f9ad287a629d5e685c41fe77a22356964f2c054183bed` |
| BQSR known sites | Known INDELs, GRCh38 | `Databases/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz` | `c0c400cfd0ca5f743b9b5e0fecd6bb05ef3f146667819944bd56919c5675d77f` |
| Exome intervals | Illumina coding targets, GRCh38 | `Databases/GRCh38/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list` | `d52f74d923cffa570ab42e3ba833aa625be67ce3ac518608d50d4013b611d934` |
| Mutect2 | Panel of normals, GRCh38, PON 5210 | `Databases/GRCh38/MuTect2.PON.5210.vcf.gz` | `2504fc52172971324c818b2e20b71cbb68ac697c6c27a2f7ecc19ec284230e35` |
| Mutect2 | gnomAD germline allele-frequency resource, GRCh38 | `Databases/GRCh38/af-only-gnomad.hg38_no_alt.vcf.gz` | `47f11c93b31f1dbc3911c528332346f570d16391de2c8b4748a08e72758a141b` |
| Structural variant calling | Delly exclusion regions, hg38 | `Databases/GRCh38/Delly/human.hg38.excl.tsv` | `caef8593f82f513694ae41b429680eab8abe39c057a119dae7760715aef006ee` |
| SAGE | Known somatic hotspots, SAGE reference bundle v5.34, GRCh38 | `SAGE/v5_34/ref/38/variants/KnownHotspots.somatic.38.vcf.gz` | `dc047d77da9e94a4f158d4f8d7b585b7cff9672a582898ed81ba193086ce9178` |
| SAGE | Actionable coding panel, SAGE reference bundle v5.34, GRCh38 | `SAGE/v5_34/ref/38/variants/ActionableCodingPanel.38.bed.gz` | `507eda10e3897f3f42973c3f0279448ff961f8e9c30781729b8c47a0f2457555` |
| SAGE | GIAB high-confidence regions v3.3.2, GRCh38 | `SAGE/v5_34/ref/38/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz` | `bd38220c2c8cadfa4f9b8ba2cfb50dc1b38fec758fb3b79cb44c16fdce59635b` |
| Reference genome | GRCh37 FASTA | `GRCh37_ref/GRCh37.fa` | `ecfa54a494fd070d7675db07ba15bdbdde4145045c2982b958d3e3553bd5bb55` |
| Genomic intervals | GRCh37 chromosome BED | `GRCh37_ref/GRCh37_chr.bed.gz` | `ad9b27f44d7e732ae3dcef026c6d236f2e5aca58629d812e126d659c059cac42` |
| Sample identity/QC | Conpair 1000 Genomes Phase 3 markers, GRCh37, 20130502, MAF 0.4, LD 0.8 | `Databases/GRCh37/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed` | `6f2dcb9668dae84731d8a66788de6c5a56258095638a29eb75fa4ee7ef00e470` |
| Sample identity/QC | Conpair 1000 Genomes Phase 3 genotypes, GRCh37, 20130502, MAF 0.4, LD 0.8 | `Databases/GRCh37/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt` | `ec7cc04d8443042f32fa11ccc1e941985d704e7d543941d76d845c16e1dfb5f3` |
| Coverage analysis | GRCh37 exome target BED | `Databases/GRCh37/GRCh37_exome.bed` | `d1238cd909b939260cdeb2e7b3ab42beac39cd0cb002475094f8098b0f6bcc71` |
| Known variants | dbSNP, assembly19 | `Databases/GRCh37/Homo_sapiens_assembly19.dbsnp.vcf` | `271eedef47bbfe228974a7b51529d9141bd8dc7cf586f71e7434715369b81bfd` |
| BQSR known sites | dbSNP, assembly19 | `Databases/GRCh37/Homo_sapiens_assembly19.dbsnp.vcf.gz` | `8ab9cd1b2632c7656ce73d78c4f2a71b2cd298a8d73acfbf18becc529ef52e2a` |
| BQSR known sites | Known INDELs, assembly19 | `Databases/GRCh37/Homo_sapiens_assembly19.known_indels.vcf.gz` | `4730fa11e61615d4056f9f1689178558979ad4f4233c9606e5f725aee3a61103` |
| Exome intervals | Broad human exome intervals, GRCh37 | `Databases/GRCh37/intervals_Broad.human.exome.b37.interval_list` | `9de6794c0d088458cd3384b9fe2cb43c9de6138a130c82f335559797a629fabc` |
| Mutect2 | WGS panel of normals, GRCh37 | `Databases/GRCh37/Mutect2-WGS-panel-b37.vcf` | `d69ace916acaab80ff2f73e6abae3a716a07d87bb9ed867b765367870c9650fe` |
| Mutect2 | WES panel of normals, GRCh37 | `Databases/GRCh37/Mutect2-exome-panel_b37.vcf` | `390b1a924328e5a6e20dd6b5df024e74cf063de926d8b753889182773b2bfc4a` |
| Mutect2 | gnomAD germline allele-frequency resource, GRCh37 | `Databases/GRCh37/af-only-gnomad.raw.sites.grch37.vcf.gz` | `b26592f0f697c6b0230c4ecf98b2307dc04f52e7d866a447bd1efe5768bd8704` |
| Structural variant calling | Delly exclusion regions, hg19 | `Databases/GRCh37/Delly/human.hg19.excl.tsv` | `249fc077cc10bc2a1481b20b32c5bf5847484ba38f9f7791bacfd51a419a7e18` |
| SAGE | Known somatic hotspots, SAGE reference bundle v5.34, GRCh37 | `SAGE/v5_34/ref/37/dna/variants/KnownHotspots.somatic.37.vcf.gz` | `cad9a9198e1d9aa9434de85caeb67f97d2899e143212238c1ec1f30ba9ba80dd` |
| SAGE | Actionable coding panel, SAGE reference bundle v5.34, GRCh37 | `SAGE/v5_34/ref/37/dna/variants/ActionableCodingPanel.37.bed.gz` | `2a6348540a2d360f6070b4713b9c54d2fcbf84fd8603f109f987c5a59343e8c9` |
| SAGE | GIAB high-confidence regions v3.2.2, GRCh37 | `SAGE/v5_34/ref/37/dna/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz` | `29257b618e908e12862513dd9324e34d883bc622d30842e55d16796965cd20f4` |
| Reference genome | mm39 FASTA | `mm39_ref/mm39.fa` | `34a1e72bb3edce582b3d5ab40fc6082039ac24fd40cf6472d73a69284c60a2b4` |
| Genomic intervals | mm39 chromosome BED | `mm39_ref/mm39_chr.bed.gz` | `9ff1050fd406d81739a5e725e9ae8de466cd67d99f9e7326dacdb9cd737bb820` |
| Coverage analysis | mm39 exome target BED | `Databases/mm39/mm39_exome.bed` | `283f7e91c111010843044997721e9bacf39f1709125ff020ddcce8418041bb21` |
| Exome intervals | mm39 exome interval list | `Databases/mm39/mm39_exome.interval_list` | `ccacf64d352cb85d53c22aa90af6823d1c72bf43d03eb54b0c4cd0869489e220` |
| Known variants | Mouse Genomes Project variant resource, mm39 | `Databases/mm39/af_only_mgp_mm39_unique.vcf.gz` | `6d97727ed06e7ef4514c69eaf4a2595eb2ed3f69b537288a7437a35066094fd6` |
| Mutect2 | Panel of normals, mm39 | `Databases/mm39/PoN.mm39.vcf.gz` | `1ceadefd7e2ae263b46aa1f7042d66b79a0a40beec8e1abb732dc9fe84d0ac8d` |
| Mutect2 | Mouse Genomes Project germline allele-frequency resource, mm39 | `Databases/mm39/af-only-mgp.mm39.vcf.gz` | `65a909b2937b97e1f227b37b54aea9bb20ea9a6c8f1f31e24756426b92850898` |
| Reference genome | rn7 FASTA | `RN7_ref/rn7.fa` | `96904dd4f449baac5f95286c67944f8104a06c9d765c4cb00380a9206564a9bd` |
| Genomic intervals | rn7 chromosome BED | `RN7_ref/rn7_chr.bed.gz` | `8f688d112e5cdd381187ad0199180ef3b495f5c2e8a558d8eb0084d736ae15d3` |
| Coverage analysis | rn7 exome target BED | `Databases/RN7/rn7_exome.bed` | `311d68190244a5f6bd9472ad8409471cf45c2198ca2e38a072b732898386399b` |
| Exome intervals | rn7 exome interval list | `Databases/RN7/rn7_exome.interval_list` | `5ddad07e8be1a60005505cd3fdda3aa0a9d52adecf02e666d5cff21f2b1980a2` |
| Mutect2 | Germline allele-frequency resource, rn7 | `Databases/RN7/af-only_rn7.vcf.gz` | `0840c6b3a620c3b0f34007b79d70b93a40b10f315d44ca5c6d854e91a3cfc1df` |

### Multi-file Database Resources

The following resources contain multiple component files. Their
individual file-level checksums are listed in
`checksums/database_files.sha256`.

| Resource | Assembly/release | Relative directory |
| -------- | ---------------- | ------------------ |
| Reference FASTA indexes and sequence dictionaries | GRCh38 | `GRCh38_ref/` |
| Reference FASTA indexes and sequence dictionaries | GRCh37 | `GRCh37_ref/` |
| Reference FASTA indexes and sequence dictionaries | mm39 | `mm39_ref/` |
| Reference FASTA indexes and sequence dictionaries | rn7 | `RN7_ref/` |
| Mutect2 20-way interval shards | GRCh38 | `Databases/GRCh38/GRCh38_interval_list_20/` |
| Mutect2 20-way interval shards | GRCh37 | `Databases/GRCh37/GRCh37_interval_list_20/` |
| Mutect2 20-way interval shards | mm39 | `Databases/mm39/mm39_interval_list_20/` |
| Mutect2 20-way interval shards | rn7 | `Databases/RN7/RN7_interval_list_20/` |
| ASCAT 1000 Genomes WGS reference files | hg38 | `Databases/GRCh38/ASCAT/WGS/hg38/` |
| ASCAT 1000 Genomes WES reference files | hg38 | `Databases/GRCh38/ASCAT/WES/hg38/` |
| ASCAT 1000 Genomes WGS reference files | hg19 | `Databases/GRCh37/ASCAT/WGS/hg19/` |
| ASCAT 1000 Genomes WES reference files | hg19 | `Databases/GRCh37/ASCAT/WES/hg19/` |
| SAGE reference bundle | v5.34, GRCh38 | `SAGE/v5_34/ref/38/` |
| SAGE reference bundle | v5.34, GRCh37 | `SAGE/v5_34/ref/37/` |

### Database Checksum Manifest

| Manifest | Files covered | SHA256 |
| -------- | ------------- | ------ |
| `checksums/database_files.sha256` | 628 | `160f66be24fa90d3d705a176c3cc1108aa0158d7b24a5d2c7850fb57d6f69f30` |

The checksum manifest was generated from a lexicographically sorted list
of database files. It records relative paths so that the database bundle
can be installed in any location.

A database installation can be verified with:

```bash
cd /path/to/database_bundle
sha256sum --check /path/to/SomaticBunny/checksums/database_files.sha256