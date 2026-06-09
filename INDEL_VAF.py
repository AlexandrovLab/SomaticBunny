import sys
import pysam


def count_indel_alleles(bam, chrom, pos0, ref, alt, min_mapq=20, min_baseq=20):
    """
    Count REF- and ALT-supporting reads for a single INDEL using CIGAR walking.

    VCF convention (standard left-anchored):
      pos0  : 0-based position of the shared anchor base (ref[0] == alt[0]).
      Deletion : REF=ACGT ALT=A    -> del_len=3, deleted bases at pos0+1..pos0+3
      Insertion: REF=A    ALT=ATCG -> ins_seq="TCG", inserted bases at pos0+1

    Returns
    -------
    (ref_count, alt_count, total, vaf, flag)
      flag  : "OK" | "MNP" | "NO_ANCHOR" | "INVALID"
      total : ref_count + alt_count
              (ambiguous / other reads are excluded from the denominator)
    """

    ref_count   = 0
    alt_count   = 0
    other_count = 0

    # ── Guard: empty alleles (malformed VCF) ─────────────────────────────────
    if not ref or not alt:
        return 0, 0, 0, 0.0, "INVALID"

    # ── Classify variant ──────────────────────────────────────────────────────
    is_del = len(ref) > len(alt)
    is_ins = len(alt) > len(ref)
    is_mnp = (not is_del) and (not is_ins) and len(ref) > 1

    if is_mnp:
        return 0, 0, 0, 0.0, "MNP"

    # Standard VCF INDELs must share the first base as the anchor
    if ref[0] != alt[0]:
        return 0, 0, 0, 0.0, "NO_ANCHOR"

    # del_len : bases deleted after the anchor   (e.g. ACGT->A  = 3)
    # ins_seq : sequence inserted after anchor   (e.g. A->ATCG  = "TCG")
    # alt[1:] is correct for standard VCF which always has exactly 1 anchor base
    del_len = len(ref) - len(alt)
    ins_seq = alt[1:] if is_ins else ""

    # ── Fetch window ──────────────────────────────────────────────────────────
    # Covers anchor + full ref/alt span + 1 bp beyond.
    # The -1 pad ensures reads whose first aligned base is exactly pos0 are fetched.
    fetch_start = max(0, pos0 - 1)
    fetch_end   = pos0 + max(len(ref), len(alt)) + 1

    for read in bam.fetch(chrom, fetch_start, fetch_end):

        # ── Read-level filters ────────────────────────────────────────────────
        if (read.is_unmapped
                or read.is_duplicate
                or read.is_secondary
                or read.is_supplementary
                or read.mapping_quality < min_mapq):
            continue

        # query_sequence can be None for certain BAM records
        if read.query_sequence is None:
            continue

        # ── Locate anchor base ────────────────────────────────────────────────
        # get_aligned_pairs(matches_only=False) yields (query_pos, ref_pos).
        # Deleted bases  → query_pos is None.
        # Inserted bases → ref_pos   is None.
        pairs = read.get_aligned_pairs(matches_only=False)

        anchor_idx  = None
        anchor_qpos = None
        for idx, (qpos, rpos) in enumerate(pairs):
            if rpos == pos0:
                anchor_idx  = idx
                anchor_qpos = qpos
                break

        # Read does not cover the anchor position at all
        if anchor_idx is None:
            continue

        # Anchor base is itself deleted in this read (pathological alignment)
        if anchor_qpos is None:
            continue

        # ── Anchor base quality filter ────────────────────────────────────────
        if read.query_qualities is not None:
            if read.query_qualities[anchor_qpos] < min_baseq:
                continue

        # ── Allele assignment ─────────────────────────────────────────────────
        if is_del:
            # ALT: exactly del_len consecutive pairs with query_pos=None after anchor
            deletion_pairs = pairs[anchor_idx + 1 : anchor_idx + 1 + del_len]

            if (len(deletion_pairs) == del_len
                    and all(qp is None for qp, rp in deletion_pairs)):
                alt_count += 1

            else:
                # Call REF only if the read spans the full deleted region
                # (del_len pairs after anchor) plus at least one base beyond,
                # confirming the read was not truncated mid-deletion.
                downstream = pairs[anchor_idx + 1 : anchor_idx + 1 + del_len + 1]
                if len(downstream) == del_len + 1:
                    ref_count += 1
                else:
                    other_count += 1   # read ends inside the deletion window

        elif is_ins:
            # ALT: exactly len(ins_seq) consecutive pairs with ref_pos=None after anchor
            ins_pairs = pairs[anchor_idx + 1 : anchor_idx + 1 + len(ins_seq)]

            if (len(ins_pairs) == len(ins_seq)
                    and all(rp is None for qp, rp in ins_pairs)):

                ins_qpositions = [qp for qp, rp in ins_pairs if qp is not None]

                if len(ins_qpositions) == len(ins_seq):
                    # Verify inserted sequence matches ALT
                    read_ins = read.query_sequence[
                        ins_qpositions[0] : ins_qpositions[-1] + 1
                    ]
                    if read_ins == ins_seq:
                        alt_count += 1
                    else:
                        other_count += 1   # insertion of same length, different seq
                else:
                    other_count += 1       # some inserted bases have no query position

            else:
                # No insertion at anchor+1 → REF, provided the read continues past anchor
                downstream = pairs[anchor_idx + 1 : anchor_idx + 2]
                if len(downstream) > 0:
                    ref_count += 1
                else:
                    other_count += 1       # read ends exactly at the anchor base

    total = ref_count + alt_count          # other_count excluded from denominator
    vaf   = alt_count / total if total > 0 else 0.0
    return ref_count, alt_count, total, vaf, "OK"


# ── Main ──────────────────────────────────────────────────────────────────────

vcf_file    = sys.argv[1]
tumor_bam   = sys.argv[2]
normal_bam  = sys.argv[3]
sample_name = sys.argv[4]
out_file    = sys.argv[5]

SKIP_FLAGS = {"MNP", "NO_ANCHOR", "INVALID"}

bam_t = pysam.AlignmentFile(tumor_bam,  "rb")
bam_n = pysam.AlignmentFile(normal_bam, "rb")

try:
    header = "\t".join([
        "sample", "#CHR", "LOC", "SAMPLE", "REF", "ALT",
        "CALLERS", "FILTER", "DOT", "INFO",
        "FORMAT", "TUMOR", "NORMAL", "VAF"
    ])

    rows = []
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            # Skip empty or malformed lines (need at least CHROM POS ID REF ALT)
            if len(fields) < 5:
                continue

            chrom   = fields[0]
            pos1    = int(fields[1])
            pos0    = pos1 - 1                  # VCF is 1-based; pysam is 0-based
            sid     = fields[2]
            ref     = fields[3]
            alt     = fields[4].split(",")[0]   # first ALT only for multi-allelic sites
            callers = fields[5] if len(fields) > 5 else "."
            filt    = fields[6] if len(fields) > 6 else "."
            dot     = fields[7] if len(fields) > 7 else "."
            info    = fields[8] if len(fields) > 8 else "."

            t_ref, t_alt, t_tot, t_vaf, t_flag = count_indel_alleles(
                bam_t, chrom, pos0, ref, alt
            )
            n_ref, n_alt, n_tot, n_vaf, n_flag = count_indel_alleles(
                bam_n, chrom, pos0, ref, alt
            )

            if t_flag in SKIP_FLAGS or n_flag in SKIP_FLAGS:
                continue

            tumor_fmt  = f"{t_ref}:{t_alt}:{t_tot}"
            normal_fmt = f"{n_ref}:{n_alt}:{n_tot}"

            rows.append("\t".join([
                sample_name, chrom, str(pos1), sid, ref, alt,
                callers, filt, dot, info,
                "RD:AD:TD", tumor_fmt, normal_fmt, f"{t_vaf:.4f}"
            ]))

    with open(out_file, "w") as out:
        out.write(header + "\n")
        if rows:
            out.write("\n".join(rows) + "\n")

finally:
    bam_t.close()
    bam_n.close()
