# EVCpipeline 
## Workflow Introduction
<img src="https://github.com/Tina04021997/EVCpipeline/blob/main/workflow_logo/v0.1.png" width="95%" height="95%">

## How to run EVC pipeline
1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment
2. Make sure you have the following folders and files in your working directory:
  - main.nf
  - nextflow.config
  - conf/base.config
  - sample.csv
  - summary.py

  > You can find these shared files and folder in `/tscc/projects/ps-lalexandrov/shared/EVC_nextflow`
3. Prepare your sample.csv file:
```
patient,status,fastq_1,fastq_2,type,file
RADS10,normal,/path/to/RADS10_normal_1.fastq.gz,/path/to/RADS10_normal_2.fastq.gz,exome,fastq
RADS10,tumor,/path/to/RADS10_tumor_1.fastq.gz,/path/to/RADS10_tumor_2.fastq.gz,exome,fastq
```
> If you are running whole exome files, specify `exome` for type and `genome` for whole genome samples in the sample.csv
> You can use the `create_sample_csv.sh` to create your sample.csv

4. Due to TSCC memory issue, you may need to modify the temp folder path to some folder in restricted:
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
nextflow run main.nf

# If your pipeline terminates with external error, or the interactive node is killed, you can resume your task after setting up the previous steps again with the following command:
nextflow run main.nf -resume

# Optionally, you can recieve an notifiction email on completion with -N flag:
nextflow run main.nf -N your_email@gmail.com
```
6. Every process result and report will be stored in the **RESULT** folder

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
