# Changelog

All notable changes to SomaticBunny are recorded here. The format follows Keep a Changelog; versions follow semantic versioning.

## [Unreleased]

## [4.0.0] - 2026-11-06

Release corresponding to the manuscript. Changes made in response to the external code audit:

- Consensus restricted to one documented region set for all callers (`--consensus_contigs`, `--callable_bed`); off-region records logged.
- SAGE frozen to 4.2 with HMF pipeline resources 2.3.0--2; jar hash recorded in every SAGE log.
- Custom reference: packaged resources validated by contig name and length, every mismatch fatal; ASCAT loci included; `--mutect2_germline` and `--mutect2_pon` overrides.
- Startup validation of `--first_step`, genome x step combinations, sample identifiers, one FASTQ pair per sample; RN7 gated behind `--allow_unvalidated_genome`.
- Results, log and temporary directories moved from the pipeline directory to the launch directory (`--outdir`).
- Whole-contig ApplyBQSR interval lists generated at runtime; complete database manifest with component table and `bin/verify_database.sh`, checked at startup.
- `--mkdup_optical_pixel_distance` (default 100).
- BSD 2-Clause license; `.gitignore`; release metadata (`CITATION.cff`, Nextflow manifest, commit in the run log).
