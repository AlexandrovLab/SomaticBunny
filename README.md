```markdown
# SMURFS: Somatic MUtation Recognition Framework & Suite

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.3-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
- [License](#license)

## Overview

The Somatic MUtation Recognition Framework & Suite (SMURFS) integrates multiple variant calling tools to improve detection accuracy for:
- Single Nucleotide Variants (SNVs)
- Insertions and Deletions (INDELs)
- Structural Variants (SVs)
- Copy Number Variants (CNVs)

This pipeline is designed to work with paired tumor-normal whole genome or whole exome sequencing data.

## Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥21.10.3
- High-performance computing environment (TSCC)
- Conda/Mamba for environment management

## Installation

1. Clone this repository or set up a working directory with the required files:

```bash
mkdir -p smurfs_pipeline/conf
cd smurfs_pipeline
```

2. Add the following required files to your working directory:
   - `main.nf`
   - `nextflow.config`
   - `conf/base.config`
   - `sample.csv`

   > **Note:** These files can be found in the shared directory: `/tscc/projects/ps-lalexandrov/shared/EVC_nextflow`

3. Prepare your `sample.csv` file according to your analysis requirements:

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
nextflow run main_stepwise_cnsv_test.nf --type [genome|exome] --step variant_calling --tool [tool1,tool2,...]
```

#### Available Tools

| Category | Available Tools |
|----------|----------------|
| SNV/INDEL Callers | `strelka2`, `mutect2`, `sage`, `muse2` |
| SV Callers | `manta`, `delly` |
| CNV Callers | `ascat`, `cnvkit` |

### Advanced Options

#### Publishing Intermediate Files

By default, intermediate files are not saved to reduce disk usage. Use the `--publish` parameter to save specific file types:

```bash
nextflow run main_conpair.nf --type exome --step variant_calling --publish raw_bam,mkdup_bam
```

**Available options:**
- `raw_bam`: Publish raw BAM files from BWA-MEM
- `recal_bam`: Publish recalibrated BAM files
- `mkdup_bam`: Publish mark duplicates BAM files

#### Resuming Failed Runs

If your pipeline terminates with an error or the interactive node is killed, resume with:

```bash
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta -resume
```

#### Email Notifications

Receive completion notification:

```bash
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta -N your_email@example.com
```

## Output

All results and reports are stored in the **RESULTS** folder with the following structure:

```
RESULTS/
├── QC/
│   ├── FastQC/
│   └── Conpair/
├── BAM/
│   ├── raw_bam/
│   ├── recal_bam/
│   └── mkdup_bam/
├── SNV_INDEL/
│   ├── Strelka2/
│   ├── Mutect2/
│   ├── SAGE/
│   └── MuSE2/
├── SV/
│   ├── Manta/
│   └── Delly/
└── CNV/
    ├── ASCAT/
    └── CNVkit/
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

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite:

```
SMURFS: A Comprehensive Somatic Variant Calling Framework (in preparation)
```
```
