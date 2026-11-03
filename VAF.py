#!/usr/bin/env python3
"""
VAF.py -- tumor / normal allele counts and VAF for every record of a SomaticBunny final consensus VCF.

Two sub-commands:

  loci   write the unique SNV loci of the VCF as a 2-column, sorted file for alleleCounter
         (alleleCounter rejects '#' lines and would count duplicated loci twice)

  table  build one fixed-schema TSV from the VCF, the two alleleCounter outputs (SNVs) and the two BAMs (indels)

SNVs   : counts come from alleleCounter, joined to the VCF on the (CHROM, POS) key.  The join is validated:
         every SNV locus must be present exactly once in each count file, and no extra locus may be present.
         REF / ALT counts are taken by base identity (Count_A/C/G/T), depth is alleleCounter's Good_depth.
INDELs : counts come from the BAMs via pysam.  The query event is first normalised (parsimonious, left-aligned
         against the reference), then a read supports ALT only if its CIGAR contains exactly that event:
         same start, same length, same inserted sequence, and the alignment continues with an aligned base
         after the event.  A read counts as REF only if it is aligned without any indel across the event
         window.  Any other read (longer / shorter / nearby indel, soft clip, read ending inside) is excluded
         from both counts.  Depth = REF + ALT.

Output columns (always 21, tab-separated, header without '#'):
  sample CHROM POS ID REF ALT TYPE FILTER CALLERS PASS_CALLERS CALLER_FILTERS
  TUMOR_REF TUMOR_ALT TUMOR_DEPTH TUMOR_VAF NORMAL_REF NORMAL_ALT NORMAL_DEPTH NORMAL_VAF STATUS NORMALIZED
FILTER        : the consensus FILTER of the VCF record (PASS | OnePassCaller | NoPassCaller)
CALLERS       : INFO callers=        (callers reporting the variant, any FILTER)
PASS_CALLERS  : INFO pass_callers=   (callers with FILTER=PASS), '.' when absent
CALLER_FILTERS: INFO caller_filters= (per-caller FILTER as caller:filter, multiple filters joined with |)
STATUS: OK | NO_COVERAGE (depth 0 in tumor or normal) | MNP | MULTIALLELIC | INVALID (not A/C/G/T)
NORMALIZED: '.' or POS:REF>ALT of the normalised indel when it differs from the VCF record.
"""
import argparse
import re
import sys

COLUMNS = ["sample", "CHROM", "POS", "ID", "REF", "ALT", "TYPE", "FILTER", "CALLERS", "PASS_CALLERS", "CALLER_FILTERS",
           "TUMOR_REF", "TUMOR_ALT", "TUMOR_DEPTH", "TUMOR_VAF",
           "NORMAL_REF", "NORMAL_ALT", "NORMAL_DEPTH", "NORMAL_VAF", "STATUS", "NORMALIZED"]
BASES = "ACGT"
NA = "NA"


def die(msg):
    sys.stderr.write("VAF.py ERROR: " + msg + "\n")
    sys.exit(1)


# ───────────────────────────── VCF ─────────────────────────────

def parse_vcf(path):
    """One dict per record of a SomaticBunny consensus VCF (sites-only, 8 columns): FILTER from column 7,
    callers / pass_callers / caller_filters from INFO. The older 9-column layout (col 6 = callers, col 7 =
    filter) and the older INFO key caller_filter= are still accepted."""
    records = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                die(f"{path} line {n}: fewer than 5 columns")
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3].upper(), f[4].upper()
            if not pos.isdigit():
                die(f"{path} line {n}: POS '{pos}' is not an integer")
            vfilter = f[6] if len(f) >= 7 and f[6] else "."
            callers, pcallers, cfilter = ".", ".", "."
            rest = "\t".join(f[5:])
            m = re.search(r"(?:^|[;\t])callers=([^;\t]*)", rest)
            if m:
                callers = m.group(1)
                m2 = re.search(r"(?:^|[;\t])pass_callers=([^;\t]*)", rest)
                pcallers = m2.group(1) if m2 else "."
                m3 = re.search(r"(?:^|[;\t])caller_filters?=([^;\t]*)", rest)
                cfilter = m3.group(1) if m3 else "."
            elif len(f) >= 7 and re.fullmatch(r"[A-Za-z,]+", f[5] or ""):
                callers, cfilter, vfilter = f[5], f[6], "."
            if "," in alt:
                vtype = "MULTIALLELIC"
            elif len(ref) == 1 and len(alt) == 1:
                vtype = "SNV"
            elif len(ref) != len(alt):
                vtype = "INDEL"
            else:
                vtype = "MNP"
            if vtype in ("SNV", "INDEL") and (set(ref) - set(BASES) or set(alt) - set(BASES)):
                vtype = "INVALID"
            records.append(dict(chrom=chrom, pos=int(pos), id=vid, ref=ref, alt=alt, type=vtype,
                                vfilter=vfilter, callers=callers, pcallers=pcallers, cfilter=cfilter, line=n))
    return records


def chrom_sort_key(c):
    s = c[3:] if c.lower().startswith("chr") else c
    return (0, int(s)) if s.isdigit() else (1, s)


# ──────────────────────── alleleCounter ────────────────────────

def read_counts(path, expected_keys):
    """alleleCounter output -> {(chrom, pos): (A, C, G, T, depth)}; validated one-to-one against expected_keys."""
    counts = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:7] != ["#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"]:
            die(f"{path}: unexpected header {header[:7]}")
        for n, line in enumerate(fh, 2):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                die(f"{path} line {n}: fewer than 7 columns")
            key = (f[0], int(f[1]))
            if key in counts:
                die(f"{path}: locus {key[0]}:{key[1]} is duplicated")
            try:
                counts[key] = tuple(int(x) for x in f[2:7])
            except ValueError:
                die(f"{path} line {n}: non-integer count")
    missing = sorted(expected_keys - counts.keys(), key=lambda k: (chrom_sort_key(k[0]), k[1]))
    extra = sorted(counts.keys() - expected_keys, key=lambda k: (chrom_sort_key(k[0]), k[1]))
    if missing:
        die(f"{path}: {len(missing)} SNV loci of the VCF are missing, first: "
            + ", ".join(f"{c}:{p}" for c, p in missing[:5]))
    if extra:
        die(f"{path}: {len(extra)} loci are not SNV loci of the VCF, first: "
            + ", ".join(f"{c}:{p}" for c, p in extra[:5]))
    return counts


def snv_counts(rec, counts):
    a, c, g, t, depth = counts[(rec["chrom"], rec["pos"])]
    by_base = dict(zip(BASES, (a, c, g, t)))
    return by_base[rec["ref"]], by_base[rec["alt"]], depth


# ───────────────────────────── indels ─────────────────────────────

def normalize_indel(fasta, chrom, pos, ref, alt):
    """Parsimonious, left-aligned representation (vt-style). pos is 1-based. Returns (pos, ref, alt)."""
    while True:
        if len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            continue
        if len(ref) == 1 and len(alt) == 1:
            break
        if ref[-1] == alt[-1] and (len(ref) == 1 or len(alt) == 1):
            # one allele would become empty: extend both to the left with the reference base
            if pos == 1:
                break
            base = fasta.fetch(chrom, pos - 2, pos - 1).upper()
            ref, alt, pos = base + ref[:-1], base + alt[:-1], pos - 1
            continue
        break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def merged_cigar(read):
    """CIGAR as a list of (op, length) with consecutive match-type ops (M, =, X) merged into 'M'."""
    out = []
    for op, ln in read.cigartuples or []:
        op = "M" if op in (0, 7, 8) else "IDNSHP"[op - 1]
        if out and out[-1][0] == "M" and op == "M":
            out[-1] = ("M", out[-1][1] + ln)
        else:
            out.append((op, ln))
    return out


def classify_read(read, pos0, del_len, ins_seq, min_baseq):
    """'ALT' | 'REF' | None for one read, for an indel whose anchor base is at 0-based pos0."""
    rpos = read.reference_start
    qpos = 0
    cig = merged_cigar(read)
    for i, (op, ln) in enumerate(cig):
        if op == "M":
            if rpos <= pos0 < rpos + ln:
                q_anchor = qpos + (pos0 - rpos)
                if read.query_qualities is not None and read.query_qualities[q_anchor] < min_baseq:
                    return None
                left_in_op = rpos + ln - 1 - pos0          # aligned bases after the anchor inside this op
                nxt = cig[i + 1] if i + 1 < len(cig) else None
                after = cig[i + 2] if i + 2 < len(cig) else None
                if del_len:                                 # deletion of del_len bases right after the anchor
                    if left_in_op == 0 and nxt and nxt[0] == "D" and nxt[1] == del_len \
                            and after and after[0] == "M":
                        return "ALT"
                    if left_in_op >= del_len + 1:            # aligned across the deleted bases and one beyond
                        return "REF"
                    return None
                else:                                       # insertion of ins_seq right after the anchor
                    if left_in_op == 0 and nxt and nxt[0] == "I" and nxt[1] == len(ins_seq) \
                            and after and after[0] == "M":
                        q_ins = qpos + ln
                        if read.query_sequence[q_ins:q_ins + len(ins_seq)].upper() == ins_seq:
                            return "ALT"
                        return None
                    if left_in_op >= 1:                      # next reference base aligned, no insertion
                        return "REF"
                    return None
            rpos += ln
            qpos += ln
        elif op == "D" or op == "N":
            if rpos <= pos0 < rpos + ln:                    # anchor itself deleted / spliced in this read
                return None
            rpos += ln
        elif op == "I" or op == "S":
            qpos += ln
        # H and P consume nothing
    return None


def indel_counts(bam, chrom, pos, ref, alt, min_mapq, min_baseq):
    if ref[0] != alt[0]:
        return None
    del_len = len(ref) - len(alt) if len(ref) > len(alt) else 0
    ins_seq = alt[1:] if len(alt) > len(ref) else ""
    pos0 = pos - 1
    n_ref = n_alt = 0
    for read in bam.fetch(chrom, max(0, pos0 - 1), pos0 + max(len(ref), len(alt)) + 1):
        if (read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary
                or read.is_qcfail or read.mapping_quality < min_mapq or read.query_sequence is None):
            continue
        cls = classify_read(read, pos0, del_len, ins_seq, min_baseq)
        if cls == "ALT":
            n_alt += 1
        elif cls == "REF":
            n_ref += 1
    return n_ref, n_alt, n_ref + n_alt


# ───────────────────────────── table ─────────────────────────────

def vaf(alt, depth):
    return f"{alt / depth:.4f}" if depth > 0 else NA


def build_table(args):
    import pysam
    records = parse_vcf(args.vcf)
    snv_keys = {(r["chrom"], r["pos"]) for r in records if r["type"] == "SNV"}
    t_counts = read_counts(args.tumor_counts, snv_keys)
    n_counts = read_counts(args.normal_counts, snv_keys)
    fasta = pysam.FastaFile(args.ref)
    bam_t = pysam.AlignmentFile(args.tumor_bam, "rb")
    bam_n = pysam.AlignmentFile(args.normal_bam, "rb")

    rows = []
    for r in records:
        na = (NA, NA, NA, NA)
        t, n, status, normalized = na, na, r["type"], "."
        if r["type"] == "SNV":
            tr, ta, td = snv_counts(r, t_counts)
            nr, nal, nd = snv_counts(r, n_counts)
            t, n = (tr, ta, td, vaf(ta, td)), (nr, nal, nd, vaf(nal, nd))
            status = "OK" if td > 0 and nd > 0 else "NO_COVERAGE"
        elif r["type"] == "INDEL":
            pos, ref, alt = normalize_indel(fasta, r["chrom"], r["pos"], r["ref"], r["alt"])
            if (pos, ref, alt) != (r["pos"], r["ref"], r["alt"]):
                normalized = f"{pos}:{ref}>{alt}"
            if ref[0] != alt[0] or len(ref) == len(alt):
                status = "INVALID"
            else:
                ct = indel_counts(bam_t, r["chrom"], pos, ref, alt, args.min_mapq, args.min_baseq)
                cn = indel_counts(bam_n, r["chrom"], pos, ref, alt, args.min_mapq, args.min_baseq)
                t, n = (ct[0], ct[1], ct[2], vaf(ct[1], ct[2])), (cn[0], cn[1], cn[2], vaf(cn[1], cn[2]))
                status = "OK" if ct[2] > 0 and cn[2] > 0 else "NO_COVERAGE"
        rows.append([args.sample, r["chrom"], str(r["pos"]), r["id"], r["ref"], r["alt"], r["type"],
                     r["vfilter"], r["callers"], r["pcallers"], r["cfilter"], *map(str, t), *map(str, n), status, normalized])

    # schema self-check before anything is written
    for row in rows:
        if len(row) != len(COLUMNS):
            die(f"internal error: row with {len(row)} columns")
        for i in (COLUMNS.index(c) for c in ("TUMOR_REF", "TUMOR_ALT", "TUMOR_DEPTH", "NORMAL_REF", "NORMAL_ALT", "NORMAL_DEPTH")):
            if row[i] != NA and not row[i].isdigit():
                die(f"internal error: non-integer count in {row[1]}:{row[2]}")
    with open(args.out, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    n_snv = sum(r["type"] == "SNV" for r in records)
    n_indel = sum(r["type"] == "INDEL" for r in records)
    sys.stderr.write(f"VAF.py: {len(rows)} records written ({n_snv} SNV, {n_indel} INDEL, "
                     f"{len(rows) - n_snv - n_indel} other) to {args.out}\n")


def write_loci(args):
    records = parse_vcf(args.vcf)
    keys = sorted({(r["chrom"], r["pos"]) for r in records if r["type"] == "SNV"},
                  key=lambda k: (chrom_sort_key(k[0]), k[1]))
    with open(args.out, "w") as out:
        for c, p in keys:
            out.write(f"{c}\t{p}\n")
    sys.stderr.write(f"VAF.py: {len(keys)} unique SNV loci written to {args.out}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("loci", help="write unique SNV loci for alleleCounter")
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("table", help="build the VAF table")
    p.add_argument("--vcf", required=True)
    p.add_argument("--tumor-counts", required=True)
    p.add_argument("--normal-counts", required=True)
    p.add_argument("--tumor-bam", required=True)
    p.add_argument("--normal-bam", required=True)
    p.add_argument("--ref", required=True, help="reference FASTA (indexed)")
    p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-mapq", type=int, default=35, help="indel reads: minimum MAPQ (alleleCounter default 35)")
    p.add_argument("--min-baseq", type=int, default=20, help="indel reads: minimum anchor base quality (alleleCounter default 20)")
    args = ap.parse_args()
    write_loci(args) if args.cmd == "loci" else build_table(args)


if __name__ == "__main__":
    main()
