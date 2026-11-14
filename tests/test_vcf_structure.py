"""
Structural checks on the consensus VCFs: every scenario file and the merged "final" VCF that
POST_SCENARIOS builds from the SNV and indel files must be parseable by htslib, declare every
INFO key they use, and be sorted.

The merge below is the same shell as in Modules/POST_scenarios.nf (header union, one #CHROM line,
records sorted by contig, position, REF, ALT); keep the two in step.
"""
import os
import subprocess

import pytest

from conftest import tool, tool_env

MERGE = r"""
set -euo pipefail
snv="$1"; indel="$2"; final="$3"
{
    awk '/^##/ && !seen[$0]++' "$snv" "$indel"
    awk '/^#CHROM/ { print; exit }' "$snv"
    awk '!/^#/' "$snv" "$indel" | sort -t$'\t' -k1,1V -k2,2n -k4,4 -k5,5
} > "$final"
"""


def merged(run, sc, tmp_path):
    final = str(tmp_path / f"scenario{sc}_final.vcf")
    subprocess.run(["bash", "-c", MERGE, "merge", run.path(sc, "snv"), run.path(sc, "indel"), final],
                   check=True, env=tool_env())
    return final


def scenario_list(run):
    return [int(s) for s in ("1,2,3,4,5,6" if run.scenarios == "all" else run.scenarios).split(",")]


def test_every_scenario_vcf_and_final_vcf_parse_with_pysam(run_default, tmp_path):
    pysam = pytest.importorskip("pysam")
    for sc in scenario_list(run_default):
        for path in (run_default.path(sc, "snv"), run_default.path(sc, "indel"), merged(run_default, sc, tmp_path)):
            vf = pysam.VariantFile(path)
            n = 0
            for rec in vf:                                    # iterating parses and validates every record
                n += 1
                assert rec.filter.keys() and set(rec.filter.keys()) <= {"PASS", "OnePassCaller", "NoPassCaller"}
                for key in rec.info.keys():
                    assert key in vf.header.info, f"{os.path.basename(path)}: INFO {key} used but not declared"
            assert n == sum(1 for l in open(path) if not l.startswith("#"))


def test_final_vcf_has_all_records_and_is_sorted(run_default, tmp_path):
    for sc in scenario_list(run_default):
        final = merged(run_default, sc, tmp_path)
        body = [l.rstrip("\n").split("\t") for l in open(final) if not l.startswith("#")]
        n_in = sum(1 for p in (run_default.path(sc, "snv"), run_default.path(sc, "indel"))
                   for l in open(p) if not l.startswith("#"))
        assert len(body) == n_in
        assert sum(1 for l in open(final) if l.startswith("#CHROM")) == 1
        keys = [(r[0], int(r[1])) for r in body]
        chrom_order = {c: i for i, c in enumerate(dict.fromkeys(r[0] for r in body))}
        assert keys == sorted(keys, key=lambda k: (chrom_order[k[0]], k[1]))


def test_bcftools_accepts_the_final_vcf(run_default, tmp_path):
    if not tool("bcftools"):
        pytest.skip("bcftools not on PATH")
    for sc in scenario_list(run_default):
        final = merged(run_default, sc, tmp_path)
        p = subprocess.run([tool("bcftools"), "view", "-H", final], capture_output=True, text=True, env=tool_env())
        assert p.returncode == 0, p.stderr
        assert "undefined tag" not in p.stderr.lower() and "not defined in the header" not in p.stderr


def contig_header(path):
    return [l.split("ID=")[1].split(",")[0] for l in open(path) if l.startswith("##contig=")]


def test_contig_header_lists_exactly_the_allowed_contigs(run_default, run_contigs):
    # default: the 24 primary chromosomes in reference order, nothing else (no chrM, viral or unplaced contig)
    assert contig_header(run_default.path(1, "snv")) == [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    # -c chr1,chr2,HPV16: only those, so the header can never name a contig a record is not allowed on
    assert contig_header(run_contigs.path(1, "snv")) == ["chr1", "chr2", "HPV16"]
