# nf-core/sarek runs

Sarek was run with `-r 3.8.0` on four datasets (WGS 30x, 60x, 90x and WES) for comparison with SomaticBunny.

## To Run

1. Edit `paths.env` (reference files, `TMPDIR`) and the `tscc_slurm` profile in `nextflow.config` (partition/account).
2. Activate a conda environment with Nextflow ≥ 25.10.2 (the benchmark used `env_nf`).
3. From an empty launch directory:

```bash
bash /path/to/benchmark/sarek/run_sarek_wgs.sh /path/to/samplesheet.csv   # WGS (same command for 30x/60x/90x)
bash /path/to/benchmark/sarek/run_sarek_wes.sh /path/to/samplesheet.csv   # WES
```

Main options: `--step mapping --aligner bwa-mem2 --genome GATK.GRCh38 --tools mutect2,strelka,freebayes,muse --only_paired_variant_calling --snv_consensus_calling --consensus_min_count 2 --normalize_vcfs`. WES adds `--wes` and `--germline_resource_tbi`, and uses the exome BED as `--intervals`.

## Configuration

Only the changes from Sarek's own configuration are included here:

- `nextflow.config`: the `tscc_slurm` profile (SLURM executor and resources), and the SNV consensus fix from Sarek 3.8.1 (`bcftools concat --allow-overlaps`, [nf-core/sarek#2128](https://github.com/nf-core/sarek/pull/2128)). The pipeline code is identical in 3.8.0 and 3.8.1, so with this fix the runs are equivalent to Sarek 3.8.1.
- `conf/base.config`: Sarek's `conf/base.config` with increased resources.

The benchmark runs used full copies of Sarek 3.8.1's `nextflow.config` and `conf/` with these changes. 

Resources in effect (resolved with `nextflow config`):

| Scope | CPUs | Memory | Sarek default |
|---|---|---|---|
| Default | 8 | 128 GB | 1 CPU, 6 GB |
| `process_low` | 16 | 250 GB | 2 CPUs, 12 GB |
| `process_medium` | 24 | 300 GB | 6 CPUs, 36 GB |
| `process_high` | 32 | 500 GB | 12 CPUs, 72 GB |

Process-specific settings are in `conf/base.config`.

## WGS intervals

`intervals/wgs_20intervals.sorted.bed` splits chr1–22, X, Y, M into 20 equal-size intervals (GATK 4.6.0.0 `SplitIntervals`), the same 20 intervals SomaticBunny uses. `intervals/make_wgs_intervals.sh` rebuilds it.

## Provenance

`pipeline_info/<run>/` holds the files Sarek wrote to `results/pipeline_info/` for each benchmark run: parameters (`params_*.json`), tool versions (`nf_core_sarek_software_mqc_versions.yml`), and execution traces, reports and timelines.

`conf/base.config` is modified from [nf-core/sarek](https://github.com/nf-core/sarek) (MIT License).