# SomaticBunny runs

SomaticBunny is run on the same four datasets as Sarek (WGS 30x, 60x, 90x and WES) for comparison with nf-core/sarek, with the scripts in this folder. Each run's records are committed under `pipeline_info/<run>/` when it finishes.

## To Run

1. Edit `paths.env` (database directory, `TMPDIR`) and, if needed, the `tscc_slurm` profile in SomaticBunny's `nextflow.config` (partition/account).
2. Activate a conda environment with Nextflow 25.10.x (the benchmark uses `env_nf`, the same environment as the Sarek runs). Nextflow 26 rejects SomaticBunny's `nextflow.config`.
3. Copy `samplesheets/samplesheet.template.csv` to `sample.csv` in an empty launch directory and fill in the FASTQ paths and sex (one tumor/normal pair per patient; ASCAT requires the `sex` column).
4. Set `SOMATICBUNNY_DIR` to the SomaticBunny checkout (the directory containing `main.nf`, at the `v1.0.0` tag) and run from that launch directory:

```bash
SOMATICBUNNY_DIR=/path/to/SomaticBunny bash /path/to/benchmark/somaticbunny/run_somaticbunny_wgs.sh   # WGS (same command for 30x/60x/90x)
SOMATICBUNNY_DIR=/path/to/SomaticBunny bash /path/to/benchmark/somaticbunny/run_somaticbunny_wes.sh   # WES
```

Extra arguments are passed through to `nextflow run` (e.g. `-resume`, or `--sample /path/to/sample.csv` for a sheet outside the launch directory).

Main options: `--first_step mapping --genome GRCh38 --publish recal_bam`, with the reference FASTA/index/dict and the whole-chromosome BED passed explicitly (see `paths.env`). WGS uses `--type genome` and `--tool ascat,manta,delly,cnvkit`. WES uses `--type exome`, `--tool ascat,manta,cnvkit` (Delly is refused for exome data by `main.nf`) and restricts the consensus to the capture targets with `--callable_bed`; the SNV/indel callers themselves are not restricted to the targets (see the `--type` note in the main README).

## Configuration

SomaticBunny is run directly from its own repository (nextflow run main.nf), so no pipeline configuration is copied here: the tscc_slurm and conda profiles (SLURM executor and resources; Conda environments with a one-hour creation timeout) are part of SomaticBunny's own nextflow.config.

The Nextflow head job is launched from an interactive SLURM session:

```bash
srun -N 1 -n 1 -c 64 --mem 900G -t 100:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash
```

with `NXF_OPTS='-Xms1g -Xmx4g'` and `TMPDIR` set to a per-run scratch directory (e.g. `.../30x` for the WGS 30x run), as encoded in the run scripts and `paths.env`.

## Reference files

All reference files come from the shared SomaticBunny database directory (`--database_path`):

| Option | File |
|---|---|
| `--ref` / `--ref_fai` / `--ref_dict` | `GRCh38_ref/GRCh38.d1.vd1.fa` (+ `.fai`, `.dict`) |
| `--bed` / `--bed_tbi` (WGS and WES) | `GRCh38_ref/hg38_chr.bed.gz` (+ `.tbi`), whole chromosomes |
| `--callable_bed` (WES) | `Databases/GRCh38/GRCh38_exome_merged.bed`, the same BED Sarek uses as `--intervals` |

## Provenance

pipeline_info/<run>/ holds, for each completed benchmark run: the exact command line (the run's .nextflow/history line, command.txt); the execution report, trace and timeline (-with-report -with-trace -with-timeline); the resolved parameters (params_<timestamp>.json, copied from RESULTS/pipeline_info/); nextflow_version.txt (written by the run scripts); the tool versions of every Conda environment the run created (conda list --explicit of each environment under work/conda/, tool_versions.txt); the database hash check (bin/verify_database.sh <database> --hash, database_check.txt); the cleanup policy (--cleanup not used, work/ kept until the disk measurement); and the disk measurement (du -sh of work/, RESULTS/ and mkdup_tmp/ at completion, disk_usage.txt). A run that has not been performed has no directory.
