# Benchmark downsampling

This directory contains the procedure used to create the downsampled tumor/normal WGS FASTQs for the SomaticBunny vs. Sarek benchmark.

## Summary

| | Input coverage | Target coverage | Fraction kept |
|---|---|---|---|
| Tumor | 90x | 60x | 0.6666 |
| Normal | 40x | 30x | 0.7500 |

Read pairs were randomly subsampled with [`seqtk sample`](https://github.com/lh3/seqtk) in fraction mode, using seed `1234` for both R1 and R2 so mates stay paired. Fractions are computed as `target / current` with `bc` at `scale=4`, which truncates rather than rounds (60/90 → 0.6666).

## Input

Paired-end tumor/normal WGS FASTQs, one set per sample, named:

```
{sample}_tumor_1.fastq.gz   {sample}_tumor_2.fastq.gz
{sample}_normal_1.fastq.gz  {sample}_normal_2.fastq.gz
```

The samples used in the benchmark are listed in [`samples.txt`](samples.txt).

## Running

Requirements: `seqtk` (benchmark run used version 1.5-r133), `gzip`, `bc`, `bash` ≥ 4.

```bash
./downsample.sh \
    -s samples.txt \
    -i /path/to/original_fastq \
    -o /path/to/fastq_downsampled_60x_or_30x \
    --jobs 4
```

## Outputs

| File | Contents |
|---|---|
| `{sample}_{tumor,normal}_{1,2}.fastq.gz` | Downsampled FASTQs |
| `run_info.txt` | Date, host, command, seqtk/gzip versions, seed, fractions |
| `manifest.tsv` | Per sample/status: fraction, seed, and R1/R2 read counts |
| `logs/{sample}.log` | Per-sample log |
