"""
Regression tests for consensus_scenarios.sh on the synthetic fixtures (tests/make_fixtures.py).

Each test names the historical failure class it guards against. Expected sets are spelled out
explicitly so that a reviewer can check them against the CASES table without running anything.
"""
import os

import pytest

C1 = "chr1"

# Variant keys (chrom, pos, ref, alt) after normalisation, by case name.
K = {
    "snv_all4_pass":      (C1, 100,  "T", "A"),
    "snv_2pass_mt_st":    (C1, 200,  "T", "A"),
    "snv_1pass_1fail":    (C1, 300,  None, None),
    "snv_muse_only":      (C1, 400,  None, None),
    "snv_all_filtered":   (C1, 500,  None, None),
    "snv_pon":            (C1, 600,  "T", "A"),
    "snv_lowtlod_rescue": (C1, 700,  "T", "A"),
    "snv_lowtlod_fail":   (C1, 800,  "G", "T"),
    "snv_dup_in_mt":      (C1, 900,  "A", "C"),
    "mnv_snv1":           (C1, 1100, "A", "G"),
    "mnv_snv2":           (C1, 1101, "C", "T"),
    "multi_alt1":         (C1, 1200, "T", "A"),
    "multi_alt2":         (C1, 1200, "T", "C"),
    "snv_offtarget_chr1": (C1, 1800, "C", "G"),
    "snv_chr2":           ("chr2", 150, "A", "C"),
    "del_equivalent":     (C1, 1000, "CA", "C"),
    "ins_2pass":          (C1, 1300, "T", "TGC"),
}


def pos_set(run, sc, kind):
    return {(r[0], r[1]) for r in run.records(sc, kind)}


def expect(names):
    return {(K[n][0], K[n][1]) for n in names}


SNV_ALL = ["snv_all4_pass", "snv_2pass_mt_st", "snv_1pass_1fail", "snv_muse_only", "snv_all_filtered", "snv_pon",
           "snv_lowtlod_rescue", "snv_lowtlod_fail", "snv_dup_in_mt", "mnv_snv1", "mnv_snv2", "multi_alt1",
           "snv_offtarget_chr1", "snv_chr2"]          # multi_alt2 shares (chrom,pos) with multi_alt1


# ---------------------------------------------------------------- scenario membership (1-vs-2 PASS support)
def test_s1_union_of_all_records_on_primary_contigs(run_default):
    assert pos_set(run_default, 1, "snv") == expect(SNV_ALL)
    assert run_default.keys(1, "indel") == {K["del_equivalent"], K["ins_2pass"]}


def test_s2_requires_one_pass_caller(run_default):
    assert pos_set(run_default, 2, "snv") == expect(SNV_ALL) - expect(["snv_all_filtered"])


def test_s3_requires_two_callers_any_filter(run_default):
    assert pos_set(run_default, 3, "snv") == expect(SNV_ALL) - expect(["snv_muse_only"])
    assert K["multi_alt2"] not in run_default.keys(3, "snv")      # second ALT: Mutect2 only


def test_s4_requires_two_callers_and_one_pass(run_default):
    assert pos_set(run_default, 4, "snv") == expect(SNV_ALL) - expect(["snv_muse_only", "snv_all_filtered"])


def test_s5_requires_two_pass_callers(run_default):
    assert pos_set(run_default, 5, "snv") == expect(SNV_ALL) - expect(["snv_muse_only", "snv_all_filtered", "snv_1pass_1fail"])
    assert run_default.decision(5, "snv", C1, 300) == "lt2PASScallers"
    assert run_default.keys(5, "indel") == {K["del_equivalent"], K["ins_2pass"]}


def test_consensus_filter_column_reflects_pass_support(run_default):
    r = {(x[0], x[1]): x for x in run_default.records(1, "snv")}
    assert r[(C1, 100)][4] == "PASS"
    assert r[(C1, 300)][4] == "OnePassCaller"
    assert r[(C1, 500)][4] == "NoPassCaller"
    assert r[(C1, 100)][5]["ncallers"] == "4" and r[(C1, 100)][5]["npass"] == "4"


# ---------------------------------------------------------------- equivalent indels
def test_equivalent_indel_representations_merge_into_one_variant(run_default):
    rec = run_default.record(5, "indel", C1, 1000)
    assert (rec[2], rec[3]) == ("CA", "C")
    assert rec[5]["ncallers"] == "3" and set(rec[5]["callers"].split(",")) == {"mt", "sa", "st"}
    # no leftover non-left-aligned spellings anywhere in scenario 1
    assert all(r[1] == 1000 for r in run_default.records(1, "indel") if r[2].startswith("CA") or r[2] == "AA")


# ---------------------------------------------------------------- duplicate records within one caller
def test_duplicate_records_in_one_caller_count_once_and_prefer_pass(run_default):
    rec = run_default.record(5, "snv", C1, 900)
    assert rec[5]["ncallers"] == "2"                       # mt counted once although it has two records
    assert "mt:PASS" in rec[5]["caller_filters"]            # the PASS record wins over the germline one
    assert rec[5].get("caller_records") == "mt:2"


# ---------------------------------------------------------------- empty caller
def test_empty_caller_file_is_handled(run_empty_muse):
    assert (C1, 400) not in pos_set(run_empty_muse, 1, "snv")          # the MuSE-only site is gone
    rec = run_empty_muse.record(5, "snv", C1, 100)
    assert rec[5]["ncallers"] == "3" and "ms" not in rec[5]["callers"].split(",")


# ---------------------------------------------------------------- MNV
def test_mnv_is_decomposed_and_the_original_is_set_aside(run_default):
    for pos in (1100, 1101):
        rec = run_default.record(5, "snv", C1, pos)
        assert set(rec[5]["callers"].split(",")) == {"sa", "st"}
        assert rec[5].get("MNV_DERIVED") is True
    mnv_file = run_default.special("sage.mnv.vcf")
    assert os.path.exists(mnv_file)
    body = [l for l in open(mnv_file) if not l.startswith("#")]
    assert len(body) == 1 and body[0].split("\t")[1:5] == ["1100", ".", "AC", "GT"]


# ---------------------------------------------------------------- multi-allelic
def test_multiallelic_record_is_split_per_alt(run_default):
    assert K["multi_alt1"] in run_default.keys(1, "snv")
    assert K["multi_alt2"] in run_default.keys(1, "snv")
    assert run_default.record(5, "snv", C1, 1200, alt="A")[5]["ncallers"] == "2"
    assert K["multi_alt2"] not in run_default.keys(5, "snv")


# ---------------------------------------------------------------- viral / non-primary contigs
def test_viral_chrM_and_unplaced_contigs_never_enter_any_scenario(run_default):
    for sc in (1, 5):
        chroms = {r[0] for r in run_default.records(sc, "snv")}
        assert chroms <= {"chr1", "chr2"}, chroms
    listing = run_default.special("outside_callable_regions.tsv")
    assert os.path.exists(listing)
    rows = [l.rstrip("\n").split("\t") for l in open(listing) if l.strip() and not l.startswith("caller")]
    excluded = {(r[0], r[2]) for r in rows}
    for contig in ("HPV16", "chrM", "chr1_KI270706v1_random"):
        assert ("mt", contig) in excluded and ("st", contig) in excluded, contig


def test_contig_allowlist_can_admit_a_viral_contig(run_contigs):
    assert ("HPV16", 150) in pos_set(run_contigs, 5, "snv")
    assert "chrM" not in {r[0] for r in run_contigs.records(1, "snv")}


# ---------------------------------------------------------------- WES target restriction
def test_callable_bed_restricts_every_caller_to_the_targets(run_targets):
    s1 = pos_set(run_targets, 1, "snv")
    assert (C1, 1800) not in s1 and ("chr2", 150) not in s1          # outside chr1:1-1500
    assert expect(["snv_all4_pass", "snv_dup_in_mt", "multi_alt1"]) <= s1
    assert run_targets.keys(1, "indel") == {K["del_equivalent"], K["ins_2pass"]}
    rows = [l.split("\t") for l in open(run_targets.special("outside_callable_regions.tsv")) if l.startswith(("mt", "st"))]
    assert any(r[2] == "chr2" for r in rows)


# ---------------------------------------------------------------- scenario 6: thresholds, rescue, PON
def test_s6_demotes_low_tlod_but_rescues_with_two_other_pass_callers(run_default, require_s6):
    rec = run_default.record(6, "snv", C1, 700)
    assert rec[4] == "PASS"
    assert "mt:lowTLOD" in rec[5]["caller_filters"] and rec[5]["npass"] == "2"


def test_s6_drops_a_variant_left_with_one_pass_caller_after_demotion(run_default, require_s6):
    assert (C1, 800) not in pos_set(run_default, 6, "snv")
    assert "lowTLOD" in run_default.decision(6, "snv", C1, 800)


def test_s6_keeps_panel_of_normals_sites_unless_veto_requested(run_default, run_pon_veto, require_s6):
    assert (C1, 600) in pos_set(run_default, 6, "snv")
    assert (C1, 600) not in pos_set(run_pon_veto, 6, "snv")
    assert "panel_of_normals" in run_pon_veto.decision(6, "snv", C1, 600)


def test_s6_bias_filter_passes_balanced_reads(run_default, require_s6):
    assert (C1, 100) in pos_set(run_default, 6, "snv")
    assert (C1, 200) in pos_set(run_default, 6, "snv")
    assert run_default.keys(6, "indel") == {K["del_equivalent"], K["ins_2pass"]}


# ---------------------------------------------------------------- argument validation
def test_invalid_scenario_and_missing_input_are_rejected(fixtures, mini_db, tmp_path):
    from conftest import run_consensus
    p = run_consensus(fixtures, mini_db, str(tmp_path / "bad_sc"), "7", expect_fail=True)
    assert p.returncode != 0 and "invalid scenario" in p.stderr
    p = run_consensus(fixtures, mini_db, str(tmp_path / "bad_muse"), "1", muse="/nonexistent.vcf", expect_fail=True)
    assert p.returncode != 0 and "not found" in p.stderr
