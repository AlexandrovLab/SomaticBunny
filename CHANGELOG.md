# Changelog

All notable changes to SomaticBunny are recorded here. The format follows Keep a Changelog; versions follow semantic versioning.

## [Unreleased]

## [1.0.1] - 2026-11-15

- CPU requests right-sized for the remaining single-threaded or lightly threaded steps (observed use in the 30x run): Delly calling and genotyping 2 CPUs / 8 GB (was 8 / 64 GB), Delly filtering 1 / 4 GB, CNVkit per tumor 8 CPUs / 48 GB (was 16 / 64 GB), ASCAT 4 / 16 GB (was 8 / 64 GB), quickchecks, report collation, renames, VAF annotation and CSV summaries 1 CPU / 4 GB (was 2 / 16 GB). No change to any tool command or output.

## [1.0.0] - 2026-11-13

Release corresponding to the manuscript. Changes made in response to the external code audit:

- Consensus restricted to one documented region set for all callers (`--consensus_contigs`, `--callable_bed`); off-region records logged.
- SAGE frozen to 4.2 with HMF pipeline resources 2.3.0--2; jar hash recorded in every SAGE log.
- Custom reference: packaged resources validated by contig name and length, every mismatch fatal; ASCAT loci included; `--mutect2_germline` and `--mutect2_pon` overrides.
- Startup validation of `--first_step`, genome x step combinations, sample identifiers, one FASTQ pair per sample; RN7 gated behind `--allow_unvalidated_genome`.
- `--tool` parsed once at startup; unknown names and tools not supported by the selected genome are refused before any job is submitted.
- mm39 and RN7 run the SomaticBunny consensus with `consensus_scenarios_nonhuman.sh` (Mutect2, Strelka2, MuSE; SNVs 2 of 3 callers, indels 2 of 2); the consensus in use is printed at startup.
- Results, log and temporary directories moved from the pipeline directory to the launch directory (`--outdir`).
- Whole-contig ApplyBQSR interval lists generated at runtime; complete database manifest with component table and `bin/verify_database.sh`, checked at startup.
- `--mkdup_optical_pixel_distance` (default 100).
- Nextflow 25.10 or later required (`manifest.nextflowVersion`).
- CI runs the consensus, VCF and VAF tests in the production `yml/POST.yml` and `yml/alleleCounter.yml` environments; callable-BED boundary test added.
- `yml/POST.yml`: unused `dkfz-bias-filter` removed, `matplotlib-base` pinned; `yml/alleleCounter.yml`: Python pinned. Stale `params.manta_bed` removed.
- Conpair documented as legacy QC; SAGE preprocessing context (no REDUX, `-skip_msi_jitter`) documented.
- BSD 2-Clause license, copyright The Regents of the University of California; `.gitignore`; release metadata (`CITATION.cff`, Nextflow manifest, commit in the run log).
- Resource requests right-sized from the peak memory and CPU use of complete runs: GetPileupSummaries 1 CPU / 8 GB (was 8 / 125 GB), BaseRecalibrator and ApplyBQSR 2 / 8 GB (was 4 / 60 GB), MarkDuplicates 2 / 48 GB (was 8 / 125 GB), bwa-mem2 64 GB (was 250 GB), MuSE 64 GB (was 250 GB), Strelka2 16 GB (was 125 GB), FastQC 2 CPUs / 8 GB (was 4 / 60 GB), Conpair 2 / 8 GB, mosdepth 4 / 8 GB, CNVkit reference 8 CPUs / 32 GB; memory and time grow with the retry attempt for every step.
- Recalibrated BAM produced directly by `samtools merge -c -p` of the 20 coordinate-sorted ApplyBQSR chunks; the second full sort (`RECALIBRATE_SortBam`) removed. One `@RG`/`@PG` per sample instead of a renamed copy per chunk. Per-sample gathers (`groupTuple(size: 20)`) release each sample as soon as its own chunks finish.
- MarkDuplicates spill files written to the task directory (node-local scratch under SLURM) instead of `<launch dir>/mkdup_tmp`; `--mkdup_temp_dir` now optional.
- Every JVM step (GATK, Picard, SAGE, Conpair, FastQC) caps its heap from the task memory and its garbage-collector threads from the task CPUs; `_JAVA_OPTIONS`/`JAVA_TOOL_OPTIONS` inherited from the submitting shell are dropped before each task.
- Startup warning when `--type exome` is used without `--callable_bed` (off-target calls retained); `--callable_bed` checked for existence. Override parameters declared before use (no "undefined parameter" warnings).
- Benchmark against nf-core/sarek moved to its own repository (`AlexandrovLab/SomaticBunnyPublication-benchmark`).