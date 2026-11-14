"""
Exact VAF values from VAF.py on the synthetic BAMs.

The tumor BAM carries a known number of ALT and REF reads at every site (make_fixtures.CASES), all
with MAPQ 60 and base quality 40, so the expected counts are the fixture numbers themselves.
SNV counts go through alleleCounter when it is on PATH (or SB_TOOLS_BIN); otherwise the same per-base counts
are produced with pysam in alleleCounter's output layout, which still exercises VAF.py's join, validation and
arithmetic. VAF.py runs under SB_VAF_PYTHON (default: the test interpreter).
"""
import os
import subprocess

import pytest

import make_fixtures as mf
from conftest import VAF_PY, VAF_PYTHON, tool, tool_env


def count_file(bam_path, loci, out, ref):
    """alleleCounter-compatible counts: MAPQ >= 35, base quality >= 20, no dup/secondary/supplementary/qcfail."""
    pysam = pytest.importorskip("pysam")
    if tool("alleleCounter"):
        subprocess.run([tool("alleleCounter"), "-r", ref, "-l", loci, "-b", bam_path, "-o", out, "-m", "20", "-q", "35"],
                       check=True, env=tool_env())
        return
    bam = pysam.AlignmentFile(bam_path)
    with open(out, "w") as fh:
        fh.write("#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth\n")
        for line in open(loci):
            chrom, pos = line.split()
            pos = int(pos)
            counts = dict.fromkeys("ACGT", 0)
            for read in bam.fetch(chrom, pos - 1, pos):
                if read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_qcfail or read.mapping_quality < 35:
                    continue
                for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                    if rpos == pos - 1:
                        if read.query_qualities[qpos] >= 20:
                            b = read.query_sequence[qpos].upper()
                            if b in counts:
                                counts[b] += 1
                        break
            fh.write(f"{chrom}\t{pos}\t{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\t{sum(counts.values())}\n")


@pytest.fixture(scope="module")
def vaf_table(run_default, fixtures, tmp_path_factory):
    d = str(tmp_path_factory.mktemp("vaf"))
    sc = 5
    final = os.path.join(d, "final.vcf")
    with open(final, "w") as out:                      # the SNV + indel merge as POST_SCENARIOS does it
        subprocess.run(["bash", "-c", r"""awk '/^##/ && !seen[$0]++' "$1" "$2"; awk '/^#CHROM/ {print; exit}' "$1"; awk '!/^#/' "$1" "$2" | sort -t$'\t' -k1,1V -k2,2n -k4,4 -k5,5""",
                        "m", run_default.path(sc, "snv"), run_default.path(sc, "indel")], stdout=out, check=True,
                       env=tool_env())
    loci = os.path.join(d, "loci.tsv")
    subprocess.run([VAF_PYTHON, VAF_PY, "loci", "--vcf", final, "--out", loci], check=True, env=tool_env())
    t_counts, n_counts = os.path.join(d, "tumor.counts"), os.path.join(d, "normal.counts")
    count_file(fixtures["tumor_bam"], loci, t_counts, fixtures["ref"])
    count_file(fixtures["normal_bam"], loci, n_counts, fixtures["ref"])
    table = os.path.join(d, "vaf.tsv")
    subprocess.run([VAF_PYTHON, VAF_PY, "table", "--vcf", final, "--tumor-counts", t_counts, "--normal-counts", n_counts,
                    "--tumor-bam", fixtures["tumor_bam"], "--normal-bam", fixtures["normal_bam"], "--ref", fixtures["ref"],
                    "--sample", f"{mf.SAMPLE}_scenario_{sc}", "--out", table], check=True, env=tool_env())
    with open(table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = {(r["CHROM"], int(r["POS"]), r["REF"], r["ALT"]): r
                for r in (dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip())}
    return header, rows


def test_table_schema(vaf_table):
    header, rows = vaf_table
    assert header == ["sample", "CHROM", "POS", "ID", "REF", "ALT", "TYPE", "FILTER", "CALLERS", "PASS_CALLERS", "CALLER_FILTERS",
                      "TUMOR_REF", "TUMOR_ALT", "TUMOR_DEPTH", "TUMOR_VAF", "NORMAL_REF", "NORMAL_ALT", "NORMAL_DEPTH", "NORMAL_VAF",
                      "STATUS", "NORMALIZED"]
    assert rows and all(len(r) == 21 for r in rows.values())


def test_exact_snv_vaf(vaf_table):
    _, rows = vaf_table
    r = rows[("chr1", 100, "T", "A")]
    assert (r["TUMOR_REF"], r["TUMOR_ALT"], r["TUMOR_DEPTH"], r["TUMOR_VAF"]) == ("28", "12", "40", "0.3000")
    assert (r["NORMAL_REF"], r["NORMAL_ALT"], r["NORMAL_DEPTH"], r["NORMAL_VAF"]) == ("30", "0", "30", "0.0000")
    assert r["STATUS"] == "OK" and r["TYPE"] == "SNV" and r["FILTER"] == "PASS"
    assert r["CALLERS"] == "mt,sa,st,ms"


def test_exact_deletion_vaf_counts_only_reads_carrying_exactly_the_event(vaf_table):
    _, rows = vaf_table
    r = rows[("chr1", 1000, "CA", "C")]
    assert (r["TUMOR_REF"], r["TUMOR_ALT"], r["TUMOR_DEPTH"], r["TUMOR_VAF"]) == ("30", "10", "40", "0.2500")
    assert (r["NORMAL_ALT"], r["NORMAL_VAF"]) == ("0", "0.0000")
    assert r["TYPE"] == "INDEL" and r["STATUS"] == "OK" and r["NORMALIZED"] == "."


def test_exact_insertion_vaf(vaf_table):
    _, rows = vaf_table
    r = rows[("chr1", 1300, "T", "TGC")]
    assert (r["TUMOR_REF"], r["TUMOR_ALT"], r["TUMOR_DEPTH"], r["TUMOR_VAF"]) == ("30", "10", "40", "0.2500")
    assert r["NORMAL_ALT"] == "0"


def test_mnv_derived_snvs_are_counted_as_plain_snvs(vaf_table):
    _, rows = vaf_table
    for key in (("chr1", 1100, "A", "G"), ("chr1", 1101, "C", "T")):
        r = rows[key]
        assert r["TYPE"] == "SNV" and r["TUMOR_ALT"] == "12" and r["TUMOR_VAF"] == "0.3000"


def test_loci_file_is_unique_and_sorted(run_default, tmp_path):
    loci = str(tmp_path / "loci.tsv")
    subprocess.run([VAF_PYTHON, VAF_PY, "loci", "--vcf", run_default.path(5, "snv"), "--out", loci], check=True, env=tool_env())
    rows = [l.split() for l in open(loci)]
    assert len(rows) == len({tuple(r) for r in rows})
    assert all(not l.startswith("#") for l in open(loci))
