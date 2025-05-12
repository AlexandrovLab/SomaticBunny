```markdown
# SMURFS: Somatic MUtation Recognition Framework & Suite

## How to run SMURFS pipeline
1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment
2. Make sure you have the following files in your working directory:
  - main.nf
  - nextflow.config
  - conf/base.config
  - sample.csv

  > You can find these files in the shared directory: `/tscc/projects/ps-lalexandrov/shared/EVC_nextflow`
3. Prepare your sample.csv file according to your analysis needs:
> For ASCAT analysis (includes sex information):
```
patient,sample,status,fastq_1,fastq_2,sex
RADS10_GBC,A,normal,/PATH/TO/NORMAL_1.fastq.gz,/PATH/TO/NORMAL_2.fastq.gz,XY
RADS10_GBC,A,tumor,/PATH/TO/TUMOR_1.fastq.gz,/PATH/TO/TUMOR_2.fastq.gz,XY
RADS10_GBC,B,tumor,/PATH/TO/TUMOR_B_1.fastq.gz,/PATH/TO/TUMOR_B_2.fastq.gz,XY
```
> Without ASCAT analysis:
```
patient,sample,status,fastq_1,fastq_2
RADS10_GBC,A,normal,/PATH/TO/NORMAL_1.fastq.gz,/PATH/TO/NORMAL_2.fastq.gz
RADS10_GBC,A,tumor,/PATH/TO/TUMOR_1.fastq.gz,/PATH/TO/TUMOR_2.fastq.gz
RADS10_GBC,B,tumor,/PATH/TO/TUMOR_B_1.fastq.gz,/PATH/TO/TUMOR_B_2.fastq.gz
```
> Specify sequencing type with `--type exome` for whole-exome or `--type genome` for whole-genome samples
4. Configure temporary directories (required for TSCC memory management):
  - Set `$params.mkdup_temp_dir` in `main.nf` (Default: `$projectDir/mkdup_tmp`)
  - Set `$workDir` in `nextflow.config` (Default: `./work`)
5. Request an interactive node and run the pipeline:

```
## Request node
srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash

## Activate Nextflow environment
conda activate env_nf

## Set TSCC temp directory
export TMPDIR=/path/to/restricted/folder/

## Run pipeline with desired tools
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta
```

## Available Tools

You can select one or more tools from each category:

- SNV/INDEL callers: strelka2, mutect2, sage, muse2
- SV callers: manta, delly
- CNV callers: ascat, cnvkit

## Publishing intermediate bam files

By default, intermediate files are not saved to reduce disk usage. To publish specific files, use the `--publish` parameter:

Available options:
- `raw_bam`: Raw BAM files from BWA_MEM
- `recal_bam`: Recalibrated BAM files 
- `mkdup_bam`: Marked duplicates BAM files

Example:
```
nextflow run main_conpair.nf --type exome --step variant_calling --publish raw_bam,mkdup_bam
```

### If your pipeline terminates with an error or the node is killed, resume with:
```
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta -resume
```

### To receive email notifications upon completion:
```
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta -N your_email@example.com
```

6. All results and reports will be stored in the **RESULTS** folder

## Tool Versions

| Tool | Version |
| --- | --- |
| FastQC | v0.12.1 |
| Picard | v2.18.27 |
| samtools | v1.21 |
| bwa-mem2 | v2.2.1 |
| Conpair | v0.2 |
| Picard MarkDuplicates | v3.2.0-1 |
| GATK4 | v4.6.0.0 |
| mosdepth | v0.3.8 |
| Strelka2 | v2.9.10 |
| Mutect2 | v4.6.0.0 (gatk) |
| SAGE | v3.3 |
| MuSE2 | v2.1.2 |
| Delly | v1.3.1 |
| CNVkit | v0.9.8 |
| ASCAT | v3.2.0 |
| Manta | v1.6.0 |
```
