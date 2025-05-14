# SMURFS: Somatic MUtation Recognition Framework & Suite

SMURFS is an ensemble pipeline for somatic variant calling that integrates multiple tools to provide comprehensive variant detection in tumor-normal paired samples.

<p align="center">
  <img src="https://github.com/AlexandrovLab/SMURFS/blob/main/workflow_logo/SMURFS.png" alt="SMURFS Pipeline Workflow"/>
</p>

## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Basic Usage](#basic-usage)
  - [Advanced Options](#advanced-options)
- [Output](#output)
- [Tool Versions](#tool-versions)

## Overview

The Somatic MUtation Recognition Framework & Suite (SMURFS) integrates multiple variant calling tools to improve detection accuracy for:
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

   > **Note:** These files can also be found in the shared directory: `/tscc/projects/ps-lalexandrov/shared/EVC_nextflow`

4. Prepare your `sample.csv` file according to your analysis requirements:

   <details>
   <summary><b>For ASCAT analysis (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,fastq_1,fastq_2,sex
   RADS10_GBC,A,normal,/PATH/TO/ERR5285401_1.fastq.gz,/PATH/TO/ERR5285401_2.fastq.gz,XY
   RADS10_GBC,A,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz,XY
   RADS10_GBC,B,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz,XY
   ```
   </details>

   <details>
   <summary><b>Without ASCAT analysis (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,fastq_1,fastq_2
   RADS10_GBC,A,normal,/PATH/TO/ERR5285401_1.fastq.gz,/PATH/TO/ERR5285401_2.fastq.gz
   RADS10_GBC,A,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz
   RADS10_GBC,B,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz
   ```
   </details>

## Configuration

For TSCC memory management, modify the following parameters:

- Set `$params.mkdup_temp_dir` in `main.nf` (Default: `$projectDir/mkdup_tmp`)
- Set `$workDir` in `nextflow.config` (Default: `./work`)

```nextflow
// In main.nf
params.mkdup_temp_dir = "/path/to/restricted/folder/mkdup_tmp"

// In nextflow.config
workDir = "/path/to/restricted/folder/work"
```

## Usage

### Basic Usage

1. Request an interactive compute node:

```bash
srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash
```

2. Activate your Nextflow environment:

```bash
conda activate env_nf
```

3. Set TSCC temporary directory:

```bash
export TMPDIR=/path/to/restricted/folder/
```

4. Run the pipeline:

```bash
nextflow run main_stepwise_cnsv_test.nf --type [genome|exome] --step [mapping|markdup|recalibration|variant_calling] --tool [ascat,manta,...]
```

#### Pipeline Steps and Tools

The SMURFS pipeline workflow is divided into multiple steps that can be run individually or in sequence using the `--step` parameter.

##### Available Steps

| Step | Description | Required Input | 
|------|-------------| ---------------|
| `mapping` | This will start from performing alignment of FASTQ files to the reference genome | fastq_1 and fastq_2 |
| `markdup` | This will start from mark duplicate reads using the aligned BAM files | raw bam and bai |
| `recalibration` | This will start from performing base quality score recalibration | mark duplicated bam and bai |
| `variant_calling` | This will start from executing variant calling tools | recalibrated bam and bai |

Example:
```bash
# Run only mapping step
nextflow run main.nf --type exome --step mapping

# Run only variant calling step
nextflow run main.nf --type exome --step variant_calling --tool strelka2,mutect2

# Run from a specific step to completion
nextflow run main.nf --type exome --step markdup --tool strelka2,mutect2,ascat
```

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
nextflow run main.nf --type exome --step variant_calling --tool manta,ascat
```

### Advanced Options

#### Publishing Intermediate Files

By default, intermediate files are not saved to reduce disk usage. Use the `--publish` parameter to save specific file types:

```bash
nextflow run main.nf --type exome --step variant_calling --publish raw_bam,mkdup_bam
```

**Available options:**
- `raw_bam`: Publish raw BAM files from BWA-MEM
- `mkdup_bam`: Publish mark duplicates BAM files
- `recal_bam`: Publish recalibrated BAM files

#### Resuming Failed Runs

If your pipeline terminates with an error or the interactive node is killed, resume with:

```bash
nextflow run main.nf --type exome --step variant_calling --tool ascat,manta -resume
```

#### Email Notifications

Receive completion notification:

```bash
nextflow run main.nf --type exome --step variant_calling --tool ascat,manta -N your_email@example.com
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
├── MuSE2
├── Mutect2
├── RECALIBRATE (optional)
├── REPORT (optional)
├── SAGE
└── STRELKA
```

## Tool Versions

| Category | Tool | Version |
|----------|------|---------|
| QC | FastQC | v0.12.1 |
| Data Processing | Picard | v2.18.27 |
| | samtools | v1.21 |
| | bwa-mem2 | v2.2.1 |
| | Conpair | v0.2 |
| | Picard MarkDuplicates | v3.2.0-1 |
| | GATK4 | v4.6.0.0 |
| | mosdepth | v0.3.8 |
| SNV/INDEL Callers | Strelka2 | v2.9.10 |
| | Mutect2 | v4.6.0.0 |
| | SAGE | v3.3 |
| | MuSE2 | v2.1.2 |
| Structural Variants | Delly | v1.3.1 |
| | Manta | v1.6.0 |
| Copy Number Variants | CNVkit | v0.9.8 |
| | ASCAT | v3.2.0 |
