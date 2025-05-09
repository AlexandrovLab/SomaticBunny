# SMURFS 
## Workflow Introduction
<img src="https://github.com/AlexandrovLab/SMURFS/blob/main/workflow_logo/v0.2.png" width="95%" height="95%">


## How to run SMURFS pipeline
1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment
2. Make sure you have the following folders and files in your working directory:
  - main.nf
  - nextflow.config
  - conf/base.config
  - sample.csv

  > You can find these shared files and folders in `/tscc/projects/ps-lalexandrov/shared/EVC_nextflow`
3. Prepare your sample.csv file:
> If you are running ASCAT:
```
patient,sample,status,fastq_1,fastq_2,sex
RADS10_GBC,A,normal,/FILE/LOCATION/ERR5285401_1.fastq.gz,/FILE/LOCATION/ERR5285401_2.fastq.gz,XY
RADS10_GBC,A,tumor,/FILE/LOCATION/ERR5285402_1.fastq.gz,/FILE/LOCATION/ERR5285402_2.fastq.gz,XY
RADS10_GBC,B,tumor,/FILE/LOCATION/ERR5285402_1.fastq.gz,/FILE/LOCATION/ERR5285402_2.fastq.gz,XY
RADS13_GBC,A,normal,/FILE/LOCATION/ERR5285404_1.fastq.gz,/FILE/LOCATION/ERR5285404_2.fastq.gz,XX
RADS13_GBC,A,tumor,/FILE/LOCATION/ERR5285405_1.fastq.gz,/FILE/LOCATION/ERR5285405_2.fastq.gz,XX
RADS13_GBC,B,tumor,/FILE/LOCATION/ERR5285405_1.fastq.gz,/FILE/LOCATION/ERR5285405_2.fastq.gz,XX
RADS17_GBC,A,normal,/FILE/LOCATION/ERR5285407_1.fastq.gz,/FILE/LOCATION/ERR5285407_2.fastq.gz,XY
RADS17_GBC,A,tumor,/FILE/LOCATION/ERR5285408_1.fastq.gz,/FILE/LOCATION/ERR5285408_2.fastq.gz,XY
```
> If you are not running ASCAT:
```
patient,sample,status,fastq_1,fastq_2
RADS10_GBC,A,normal,/FILE/LOCATION/ERR5285401_1.fastq.gz,/FILE/LOCATION/ERR5285401_2.fastq.gz
RADS10_GBC,A,tumor,/FILE/LOCATION/ERR5285402_1.fastq.gz,/FILE/LOCATION/ERR5285402_2.fastq.gz
RADS10_GBC,B,tumor,/FILE/LOCATION/ERR5285402_1.fastq.gz,/FILE/LOCATION/ERR5285402_2.fastq.gz
RADS13_GBC,A,normal,/FILE/LOCATION/ERR5285404_1.fastq.gz,/FILE/LOCATION/ERR5285404_2.fastq.gz
RADS13_GBC,A,tumor,/FILE/LOCATION/ERR5285405_1.fastq.gz,/FILE/LOCATION/ERR5285405_2.fastq.gz
RADS13_GBC,B,tumor,/FILE/LOCATION/ERR5285405_1.fastq.gz,/FILE/LOCATION/ERR5285405_2.fastq.gz
RADS17_GBC,A,normal,/FILE/LOCATION/ERR5285407_1.fastq.gz,/FILE/LOCATION/ERR5285407_2.fastq.gz
RADS17_GBC,A,tumor,/FILE/LOCATION/ERR5285408_1.fastq.gz,/FILE/LOCATION/ERR5285408_2.fastq.gz
```
> If you are running whole-exome files, specify `--type exome` and `--type genome` for whole-genome samples when running the command
```
4. Due to TSCC memory issue, you may need to modify the temp folder path to a folder in the restricted:
  - `$params.mkdup_temp_dir` in the `main.nf` file (Default: `$projectDir/mkdup_tmp`)
  -  `$workDir` in the `nextflow.config` file (Default: `./work`)
5. Request an interactive node and run Nextflow in your working directory under an interactive node:

```
# Node requesting
srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash

# Activate your nextflow conda environment
conda activate env_nf

# Export TSCC temp directory to any folder in restricted
export TMPDIR=/some/folder/in/restricted/

# Run nextflow
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta

# If your pipeline terminates with an external error, or the interactive node is killed, you can resume your task after setting up the previous steps again with the following command:
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta-resume

# Optionally, you can receive an notification email on completion with -N flag:
nextflow run main_stepwise_cnsv_test.nf --type exome --step variant_calling --tool ascat,manta -N your_email@gmail.com
```
6. Every result and report will be stored in the **RESULTS** folder

## Tool Versions

| Tool | Version |
| --- | --- |
| FastQC | v0.12.1 |
| Picard | v2.18.27 |
| samtool | v1.21 |
| bwa-mem2 | v2.2.1 |
| Conpair | v0.2 |
| Picard MarkDuplicates | v3.2.0-1 |
| gatk4 | v4.6.0.0 |
| mosdepth | v0.3.8 |
| Strelka2 | v2.9.10 |
| Mutect2 | v4.6.0.0 (gatk) |
| SAGE | v3.3 |
| MuSE2 | v2.1.2 |
| Delly | v1.3.1 |
| CNVkit | v0.9.8 |
| ASCAT | v3.2.0 |
| Manta | v1.6.0 |
