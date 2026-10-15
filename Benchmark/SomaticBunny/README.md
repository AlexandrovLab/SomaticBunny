# SomaticBunny runs

SomaticBunny was run on the same four datasets as Sarek (WGS 30x, 60x, 90x and WES) for comparison with nf-core/sarek.

## To Run

1. Edit `paths.env` (reference files, `TMPDIR`) and, if needed, the `tscc_slurm` profile in SomaticBunny's `nextflow.config` (partition/account).
2. Activate a conda environment with Nextflow ≥ 25.10.2 (the benchmark used `env_nf`, the same environment as the Sarek runs).
3. Set `SOMATICBUNNY_DIR` to the SomaticBunny checkout (the directory containing `main.nf`) and run from an empty launch directory:

```bash
SOMATICBUNNY_DIR=/path/to/SomaticBunny bash /path/to/benchmark/somaticbunny/run_somaticbunny_wgs.sh   # WGS (same command for 30x/60x/90x)
SOMATICBUNNY_DIR=/path/to/SomaticBunny bash /path/to/benchmark/somaticbunny/run_somaticbunny_wes.sh   # WES
```

Main options: `--first_step mapping --genome GRCh38 --publish recal_bam --tool ascat,manta,delly,cnvkit`, with the reference FASTA/index/dict and genome BED passed explicitly (see `paths.env`). WGS uses `--type genome`; WES uses `--type exome` and the exome BED.

## Configuration

SomaticBunny was run directly from its own repository (`nextflow run main.nf`), so no pipeline configuration is copied here — the `tscc_slurm` profile (SLURM executor and resources) is part of SomaticBunny's own `nextflow.config`.

The Nextflow head job was launched from an interactive SLURM session:

```bash
srun -N 1 -n 1 -c 64 --mem 900G -t 100:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash
```

with `NXF_OPTS='-Xms1g -Xmx4g'` and `TMPDIR` set to a per-run scratch directory (e.g. `.../30x` for the WGS 30x run), as encoded in the run scripts and `paths.env`.

## Reference files

All reference files come from the shared `EVC_nextflow` database directory (`--database_path`):

| Option | File |
|---|---|
| `--ref` / `--ref_fai` / `--ref_dict` | `GRCh38_ref/GRCh38.d1.vd1.fa` (+ `.fai`, `.dict`) |
| `--bed` / `--bed_tbi` (WGS) | `GRCh38_ref/hg38_chr.bed.gz` (+ `.tbi`) |
| `--bed` / `--bed_tbi` (WES) | exome BED (see `paths.env`) |

## Provenance

`pipeline_info/<run>/` holds the files SomaticBunny/Nextflow wrote for each benchmark run: execution reports (`-with-report`), traces and timelines, and the exact command lines. `nextflow_version.txt` (written by the run scripts) records the Nextflow version used.
