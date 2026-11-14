# Changelog

All notable changes to SomaticBunny are recorded here. The format follows Keep a Changelog; versions follow semantic versioning.

## [Unreleased]

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
- Benchmark scripts replaced with the v1.0.0 protocol.