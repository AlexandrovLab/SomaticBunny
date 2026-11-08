#!/usr/bin/env python3
"""
Synthetic, deterministic test data for SomaticBunny's consensus and VAF steps.

Everything is generated from a seeded random reference, so no binary fixture is committed:
  ref/test.fa (+ .fai, .dict)   all 24 human primary chromosomes (tiny), chrM, a viral contig (HPV16)
                                and an unplaced contig (chr1_KI270706v1_random)
  vcf/<sample>.<caller>.vcf[.gz] one VCF per caller in that caller's own layout (sample order, INFO keys)
  bam/tumor.bam, bam/normal.bam reads with known REF/ALT counts at every variant site
  targets.bed                   a "capture kit" covering chr1:1-1500 only
  bseq_header                   header line prepended by consensus_scenarios.sh before the bias filter

The CASES table below is the single source of truth: tests import it to derive their expectations.
usage: make_fixtures.py OUT_DIR
"""
import os
import random
import subprocess
import sys

SAMPLE = "T001_A"
READ_LEN = 50
CALLERS = ("mt", "ms", "st", "sa")          # Mutect2, MuSE, Strelka2, SAGE (consensus tags)

# ------------------------------------------------------------------ reference
CONTIGS = [("chr1", 3000)] + [(f"chr{i}", 300) for i in range(2, 23)] + [("chrX", 300), ("chrY", 300),
           ("chrM", 300), ("HPV16", 400), ("chr1_KI270706v1_random", 300)]


def make_reference():
    rng = random.Random(20261106)
    seqs = {}
    for name, length in CONTIGS:
        seqs[name] = [rng.choice("ACGT") for _ in range(length)]
    c1 = seqs["chr1"]
    # homopolymer for the equivalent-indel case: 999=G 1000=C 1001-1005=A 1006=T  (1-based)
    for pos, base in zip(range(999, 1008), "GCAAAAATG"):
        c1[pos - 1] = base
    # unambiguous insertion anchor for the insertion case: 1299=A 1300=T 1301=A 1302=C
    for pos, base in zip(range(1298, 1304), "CATACG"):
        c1[pos - 1] = base
    # MNV case needs two bases at 1100-1101 that differ from their ALT: fix them to A C
    c1[1100 - 1], c1[1101 - 1] = "A", "C"
    return {k: "".join(v) for k, v in seqs.items()}


def alt_of(base, n=1):
    """A deterministic ALT base: the n-th next base in ACGT order."""
    return "ACGT"[("ACGT".index(base) + n) % 4]


# ------------------------------------------------------------------ cases
# Each case: name, contig, pos, kind, per-caller records {tag: (FILTER, score, extra)} where score is
# TLOD (mt) or SomaticEVS (st); None = caller did not report it. 'ref'/'alt' are filled from the sequence
# for SNVs; indels and the MNV give them explicitly. 'reads' = (n_alt, n_ref) in the tumor BAM.
CASES = [
    dict(name="snv_all4_pass",      contig="chr1", pos=100,  kind="snv",
         calls=dict(mt=("PASS", 35.0), ms=("PASS", None), st=("PASS", 25.0), sa=("PASS", None)), reads=(12, 28)),
    dict(name="snv_2pass_mt_st",    contig="chr1", pos=200,  kind="snv",
         calls=dict(mt=("PASS", 30.0), st=("PASS", 20.0)), reads=(12, 28)),
    dict(name="snv_1pass_1fail",    contig="chr1", pos=300,  kind="snv",
         calls=dict(mt=("PASS", 30.0), sa=("minTumorQual", None)), reads=(12, 28)),
    dict(name="snv_muse_only",      contig="chr1", pos=400,  kind="snv",
         calls=dict(ms=("PASS", None)), reads=(12, 28)),
    dict(name="snv_all_filtered",   contig="chr1", pos=500,  kind="snv",
         calls=dict(mt=("weak_evidence", 4.0), st=("LowEVS", 5.0)), reads=(3, 37)),
    dict(name="snv_pon",            contig="chr1", pos=600,  kind="snv",
         calls=dict(mt=("panel_of_normals", 40.0), st=("PASS", 25.0), sa=("PASS", None)), reads=(12, 28)),
    dict(name="snv_lowtlod_rescue", contig="chr1", pos=700,  kind="snv",
         calls=dict(mt=("PASS", 5.0), st=("PASS", 20.0), sa=("PASS", None)), reads=(12, 28)),
    dict(name="snv_lowtlod_fail",   contig="chr1", pos=800,  kind="snv",
         calls=dict(mt=("PASS", 5.0), st=("PASS", 20.0)), reads=(12, 28)),
    dict(name="snv_dup_in_mt",      contig="chr1", pos=900,  kind="snv",
         calls=dict(mt=[("germline", 8.0), ("PASS", 25.0)], st=("PASS", 22.0)), reads=(12, 28)),
    dict(name="mnv_sage",           contig="chr1", pos=1100, kind="mnv", ref="AC", alt="GT",
         calls=dict(sa=("PASS", None)), reads=(12, 28)),            # Strelka calls the two SNVs separately
    dict(name="mnv_snv1_strelka",   contig="chr1", pos=1100, kind="snv", ref="A", alt="G",
         calls=dict(st=("PASS", 21.0)), reads=None),
    dict(name="mnv_snv2_strelka",   contig="chr1", pos=1101, kind="snv", ref="C", alt="T",
         calls=dict(st=("PASS", 21.0)), reads=None),
    dict(name="multiallelic_mt",    contig="chr1", pos=1200, kind="multi",
         calls=dict(mt=("PASS", (33.0, 12.0)), ms=("PASS", None)), reads=(12, 28)),   # MuSE reports ALT1 only
    dict(name="del_equivalent",     contig="chr1", pos=1000, kind="del",            # CA>C at 1000 (left-aligned)
         calls=dict(mt=("PASS", 28.0), st=("PASS", 24.0), sa=("PASS", None)), reads=(10, 30)),
    dict(name="ins_2pass",          contig="chr1", pos=1300, kind="ins", ref="T", alt="TGC",
         calls=dict(mt=("PASS", 26.0), sa=("PASS", None)), reads=(10, 30)),
    dict(name="snv_offtarget_chr1", contig="chr1", pos=1800, kind="snv",
         calls=dict(mt=("PASS", 30.0), st=("PASS", 20.0)), reads=(12, 28)),
    dict(name="snv_chr2",           contig="chr2", pos=150,  kind="snv",
         calls=dict(mt=("PASS", 30.0), st=("PASS", 20.0)), reads=(12, 28)),
    dict(name="snv_viral",          contig="HPV16", pos=150, kind="snv",
         calls=dict(mt=("PASS", 30.0), st=("PASS", 20.0)), reads=(12, 28)),
    dict(name="snv_chrM",           contig="chrM", pos=50,   kind="snv",
         calls=dict(mt=("PASS", 30.0), st=("PASS", 20.0)), reads=None),
    dict(name="snv_random_contig",  contig="chr1_KI270706v1_random", pos=100, kind="snv",
         calls=dict(mt=("PASS", 30.0), st=("PASS", 20.0)), reads=None),
]

# Representations each caller uses for the equivalent deletion (all normalise to chr1:1000 CA>C).
DEL_REPR = {"mt": (1000, "CA", "C"), "st": (1004, "AA", "A"), "sa": (1001, "AA", "A")}


def resolve(seqs):
    """Fill ref/alt for SNVs and the multiallelic case from the sequence."""
    for c in CASES:
        seq = seqs[c["contig"]]
        if c["kind"] == "snv" and "ref" not in c:
            c["ref"] = seq[c["pos"] - 1]
            c["alt"] = alt_of(c["ref"])
        elif c["kind"] == "multi":
            c["ref"] = seq[c["pos"] - 1]
            c["alt"] = alt_of(c["ref"], 1) + "," + alt_of(c["ref"], 2)
        elif c["kind"] == "del":
            c["ref"], c["alt"] = "CA", "C"
        assert c["kind"] == "mnv" or c["kind"] in ("del", "ins") or c["ref"] == seq[c["pos"] - 1] or c["kind"] == "multi"
    return CASES


# ------------------------------------------------------------------ VCF writers
def contig_lines():
    return [f"##contig=<ID={n},length={l}>" for n, l in CONTIGS]


def records_for(tag):
    """Yield (contig, pos, ref, alt, FILTER, score) for one caller, in the caller's own representation."""
    out = []
    for c in CASES:
        calls = c["calls"].get(tag)
        if calls is None:
            continue
        calls = calls if isinstance(calls, list) else [calls]
        for filt, score in calls:
            if c["kind"] == "del":
                pos, ref, alt = DEL_REPR[tag]
            else:
                pos, ref, alt = c["pos"], c["ref"], c["alt"]
            out.append((c["contig"], pos, ref, alt, filt, score))
    order = {n: i for i, (n, _) in enumerate(CONTIGS)}
    out.sort(key=lambda r: (order[r[0]], r[1], r[3]))
    return out


def write_vcf(path, header_lines, columns, rows, gz):
    text = "\n".join(["##fileformat=VCFv4.2"] + header_lines + contig_lines()
                     + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(columns)]
                     + ["\t".join(map(str, r)) for r in rows]) + "\n"
    if gz:
        plain = path[:-3]
        open(plain, "w").write(text)
        subprocess.run(["bgzip", "-f", plain], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", path], check=True)
    else:
        open(path, "w").write(text)


def n_alts(alt):
    return alt.count(",") + 1


def mutect2_vcf(path):
    hdr = ["##tumor_sample=TUMOR", "##normal_sample=NORMAL",
           "##FILTER=<ID=PASS,Description=\"All filters passed\">",
           "##FILTER=<ID=germline,Description=\"Evidence indicates the site is germline\">",
           "##FILTER=<ID=weak_evidence,Description=\"Mutation does not meet likelihood threshold\">",
           "##FILTER=<ID=panel_of_normals,Description=\"Blacklisted site in panel of normals\">",
           "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">",
           "##INFO=<ID=TLOD,Number=A,Type=Float,Description=\"Log 10 likelihood ratio score of variant existing versus not existing\">",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
           "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions\">",
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"]
    rows = []
    for contig, pos, ref, alt, filt, score in records_for("mt"):
        k = n_alts(alt)
        tlod = ",".join(f"{s:.2f}" for s in (score if isinstance(score, tuple) else (score,)))
        ad_t = "28," + ",".join(["12"] * k)
        af_t = ",".join(["0.300"] * k)
        ad_n = "30," + ",".join(["0"] * k)
        af_n = ",".join(["0.000"] * k)
        gt_t = "0/1" if k == 1 else "0/1/2"
        rows.append((contig, pos, ".", ref, alt, ".", filt, f"DP=70;TLOD={tlod}", "GT:AD:AF:DP",
                     f"{gt_t}:{ad_t}:{af_t}:40", f"0/0:{ad_n}:{af_n}:30"))
    write_vcf(path, hdr, ["TUMOR", "NORMAL"], rows, gz=False)     # tumor first: exercises the column swap


def muse_vcf(path):
    hdr = ["##FILTER=<ID=PASS,Description=\"Accept as a confident somatic mutation\">",
           "##FILTER=<ID=Tier2,Description=\"Less confident somatic mutation\">",
           "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
           "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
           "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Somatic status\">"]
    rows = []
    for contig, pos, ref, alt, filt, score in records_for("ms"):
        alt1 = alt.split(",")[0]                    # MuSE reports one ALT
        rows.append((contig, pos, "rs999", ref, alt1, ".", filt, "SOMATIC", "GT:DP:AD:SS",
                     "0/1:40:28,12:2", "0/0:30:30,0:0"))
    write_vcf(path, hdr, ["TUMOR", "NORMAL"], rows, gz=False)     # MuSE writes the tumor first


def strelka_snv_vcf(path):
    hdr = ["##FILTER=<ID=PASS,Description=\"All filters passed\">",
           "##FILTER=<ID=LowEVS,Description=\"Somatic Empirical Variant Score (SomaticEVS) is below threshold\">",
           "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">",
           "##INFO=<ID=QSS,Number=1,Type=Integer,Description=\"Quality score for any somatic snv\">",
           "##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description=\"Somatic Empirical Variant Score\">",
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
           "##FORMAT=<ID=AU,Number=2,Type=Integer,Description=\"Number of 'A' alleles (tier1,tier2)\">"]
    rows = []
    for contig, pos, ref, alt, filt, score in records_for("st"):
        if len(ref) != len(alt) or len(ref) != 1:
            continue
        rows.append((contig, pos, ".", ref, alt, ".", filt, f"SOMATIC;QSS=60;SomaticEVS={score:.2f}", "DP:AU",
                     "30:0,0", "40:0,0"))
    write_vcf(path, hdr, ["NORMAL", "TUMOR"], rows, gz=True)


def strelka_indel_vcf(path):
    hdr = ["##FILTER=<ID=PASS,Description=\"All filters passed\">",
           "##FILTER=<ID=LowEVS,Description=\"Somatic Empirical Variant Score (SomaticEVS) is below threshold\">",
           "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">",
           "##INFO=<ID=QSI,Number=1,Type=Integer,Description=\"Quality score for any somatic variant\">",
           "##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description=\"Somatic Empirical Variant Score\">",
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
           "##FORMAT=<ID=TAR,Number=2,Type=Integer,Description=\"Reads strongly supporting alternate allele\">"]
    rows = []
    for contig, pos, ref, alt, filt, score in records_for("st"):
        if len(ref) == len(alt):
            continue
        rows.append((contig, pos, ".", ref, alt, ".", filt, f"SOMATIC;QSI=55;SomaticEVS={score:.2f}", "DP:TAR",
                     "30:0,0", "40:10,10"))
    write_vcf(path, hdr, ["NORMAL", "TUMOR"], rows, gz=True)


def sage_vcf(path):
    hdr = ["##FILTER=<ID=PASS,Description=\"All filters passed\">",
           "##FILTER=<ID=minTumorQual,Description=\"Insufficient tumor quality\">",
           "##INFO=<ID=TIER,Number=1,Type=String,Description=\"Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]\">",
           "##INFO=<ID=MH,Number=1,Type=String,Description=\"Microhomology\">",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"]
    rows = []
    for contig, pos, ref, alt, filt, score in records_for("sa"):
        rows.append((contig, pos, ".", ref, alt, "150", filt, "TIER=HIGH_CONFIDENCE", "GT:AD:DP",
                     "0/0:30,0:30", "0/1:28,12:40"))
    write_vcf(path, hdr, ["NORMAL", "TUMOR"], rows, gz=True)


# ------------------------------------------------------------------ BAMs
def make_bams(seqs, out_dir):
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in CONTIGS],
              "RG": [{"ID": "T001_A_tumor", "SM": "T001_A_tumor", "PL": "ILLUMINA"}]}
    tid = {n: i for i, (n, _) in enumerate(CONTIGS)}

    def read(name, contig, start0, cigar, seq, flag):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.reference_id = tid[contig]
        a.reference_start = start0
        a.mapping_quality = 60
        a.cigarstring = cigar
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.flag = flag
        a.next_reference_id = tid[contig]
        a.next_reference_start = start0
        a.template_length = 0
        a.set_tag("RG", "T001_A_tumor")
        return a

    def site_reads(c, n_alt, n_ref, prefix):
        """Reads over one site; ALT and REF reads alternate strand and read1/read2 so no bias filter fires."""
        seq = seqs[c["contig"]]
        pos0 = c["pos"] - 1
        out = []
        for i in range(n_alt + n_ref):
            is_alt = i < n_alt
            offset = 10 + (i * 7) % 25               # site sits 10..34 bases into the read
            s0 = pos0 - offset
            flag = 1 | 2 | (16 if i % 2 else 0) | (64 if (i // 2) % 2 == 0 else 128)
            if not is_alt or c["kind"] in ("snv", "multi", "mnv"):
                bases = list(seq[s0:s0 + READ_LEN])
                if is_alt:
                    for j, b in enumerate(c["alt"].split(",")[0]):     # one base, or both bases of the MNV
                        bases[offset + j] = b
                out.append(read(f"{prefix}{c['name']}_{i}", c["contig"], s0, f"{READ_LEN}M", "".join(bases), flag))
            elif c["kind"] == "del":                 # CA>C: anchor at pos0, the base after it deleted
                left = offset + 1
                cigar = f"{left}M1D{READ_LEN - left}M"
                bases = seq[s0:s0 + left] + seq[s0 + left + 1:s0 + left + 1 + READ_LEN - left]
                out.append(read(f"{prefix}{c['name']}_{i}", c["contig"], s0, cigar, bases, flag))
            elif c["kind"] == "ins":                 # T>TGC: insertion right after the anchor
                left = offset + 1
                ins = c["alt"][1:]
                cigar = f"{left}M{len(ins)}I{READ_LEN - left}M"
                bases = seq[s0:s0 + left] + ins + seq[s0 + left:s0 + READ_LEN]
                out.append(read(f"{prefix}{c['name']}_{i}", c["contig"], s0, cigar, bases, flag))
        return out

    for label, as_tumor in (("tumor", True), ("normal", False)):
        reads = []
        for c in CASES:
            if c["reads"] is None:
                continue
            n_alt, n_ref = c["reads"] if as_tumor else (0, 30)
            reads.extend(site_reads(c, n_alt, n_ref, f"{label}_"))
        reads.sort(key=lambda a: (a.reference_id, a.reference_start))
        path = os.path.join(out_dir, f"{label}.bam")
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for a in reads:
                bam.write(a)
        pysam.index(path)


# ------------------------------------------------------------------ main
def build(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for sub in ("ref", "vcf", "bam"):
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
    seqs = make_reference()
    resolve(seqs)
    fa = os.path.join(out_dir, "ref", "test.fa")
    with open(fa, "w") as fh:
        for name, _ in CONTIGS:
            fh.write(f">{name}\n")
            s = seqs[name]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    subprocess.run(["samtools", "faidx", fa], check=True)
    subprocess.run(["samtools", "dict", "-o", fa[:-3] + ".dict", fa], check=True)
    v = os.path.join(out_dir, "vcf")
    mutect2_vcf(os.path.join(v, f"{SAMPLE}.mutect2.vcf"))
    muse_vcf(os.path.join(v, f"{SAMPLE}.muse.vcf"))
    strelka_snv_vcf(os.path.join(v, f"{SAMPLE}.strelka_snv.vcf.gz"))
    strelka_indel_vcf(os.path.join(v, f"{SAMPLE}.strelka_indel.vcf.gz"))
    sage_vcf(os.path.join(v, f"{SAMPLE}.sage.vcf.gz"))
    # an "empty caller": header only
    write_vcf(os.path.join(v, "empty.muse.vcf"), ["##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"x\">"],
              ["TUMOR", "NORMAL"], [], gz=False)
    make_bams(seqs, os.path.join(out_dir, "bam"))
    with open(os.path.join(out_dir, "targets.bed"), "w") as fh:
        fh.write("chr1\t0\t1500\n")
    with open(os.path.join(out_dir, "bseq_header"), "w") as fh:
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
    return seqs


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixtures_out")
    print("fixtures written")
