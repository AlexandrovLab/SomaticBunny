#!/usr/bin/env bash
#######################################################################
# SomaticBunny benchmark: build the 20-interval BED passed to Sarek with
# --intervals for the WGS runs (wgs_20intervals.sorted.bed).
#
# chr1-22, X, Y, M are split into 20 intervals of equal total length
# with GATK SplitIntervals. These are the same 20 interval lists that
# SomaticBunny uses to scatter BQSR and Mutect2
# (Databases/GRCh38/GRCh38_interval_list_20), so both pipelines get
# the same 20-way split of the genome.
#
# Environment variables:
#   REF           Reference FASTA with .fai and .dict (benchmark: GRCh38.fa)
#   GATK          gatk launcher (default: gatk on PATH; benchmark: GATK 4.6.0.0)
#   INTERVAL_DIR  Convert existing *-scattered.interval_list files from this
#                 directory instead of running SplitIntervals, e.g. the
#                 SomaticBunny database dir Databases/GRCh38/GRCh38_interval_list_20
#   PER_INTERVAL_BED=1   Also write one BED per interval (optional, for inspection)
#   COMPARE       Existing wgs_20intervals.sorted.bed to check the output against
#######################################################################
set -euo pipefail

GATK="${GATK:-gatk}"
SCATTER=20
OUT_BED="wgs_20intervals.sorted.bed"

if [[ -n "${COMPARE:-}" ]]; then
    [[ -f "$COMPARE" ]] || { echo "ERROR: COMPARE file not found: $COMPARE" >&2; exit 1; }
    if [[ -e "$OUT_BED" && "$(realpath "$COMPARE")" == "$(realpath "$OUT_BED")" ]]; then
        OUT_BED="wgs_20intervals.regenerated.bed"
        echo "COMPARE is ${COMPARE}; writing the new BED to ${OUT_BED} instead"
    fi
fi

if [[ -z "${INTERVAL_DIR:-}" ]]; then
    : "${REF:?set REF=/path/to/GRCh38.fa (or INTERVAL_DIR to convert existing interval lists)}"
    command -v "$GATK" >/dev/null 2>&1 || { echo "ERROR: gatk not found (set GATK=/path/to/gatk)" >&2; exit 1; }

    # 1. Keep only chr1-22 and chrX, chrY, chrM
    for i in {1..22} X Y M; do echo "chr$i"; done > main_chromosomes.list

    # 2. Split into 20 intervals with equal numbers of bases
    "$GATK" SplitIntervals \
        -R "$REF" \
        -L ./main_chromosomes.list \
        --scatter-count "$SCATTER" \
        -O ./interval_list

    # 3. Rename 0000-0019-scattered.interval_list to 1-20-scattered.interval_list
    for i in $(seq 0 $((SCATTER - 1))); do
        mv "interval_list/$(printf '%04d' "$i")-scattered.interval_list" \
           "interval_list/$((i + 1))-scattered.interval_list"
    done
    INTERVAL_DIR="./interval_list"
fi

for n in $(seq 1 "$SCATTER"); do
    [[ -f "${INTERVAL_DIR}/${n}-scattered.interval_list" ]] || {
        echo "ERROR: missing ${INTERVAL_DIR}/${n}-scattered.interval_list" >&2; exit 1; }
done

# (Optional) one BED per interval. how I inspect the sections
if [[ "${PER_INTERVAL_BED:-0}" == "1" ]]; then
    for n in $(seq 1 "$SCATTER"); do
        "$GATK" IntervalListToBed -I "${INTERVAL_DIR}/${n}-scattered.interval_list" -O "${n}.bed"
    done
fi

# 4. Convert all interval lists to one BED (interval_list is 1-based,
#    BED is 0-based, hence $2-1). Files are read in order 1..20, which is
#    genome order (chr1..chr22, chrX, chrY, chrM), so the output is sorted.
for n in $(seq 1 "$SCATTER"); do
    f="${INTERVAL_DIR}/${n}-scattered.interval_list"
    grep -q -v '^@' "$f" || { echo "ERROR: no intervals in $f" >&2; exit 1; }
    grep -v '^@' "$f" | awk '{print $1 "\t" $2-1 "\t" $3}'
done > "$OUT_BED"

echo "Wrote ${OUT_BED}: $(wc -l < "$OUT_BED") regions, $(awk '{s += $3 - $2} END {print s}' "$OUT_BED") bp"

# Compressed + indexed copy (as in the original procedure; Sarek was given the
# uncompressed BED)
if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
    bgzip -c "$OUT_BED" > "${OUT_BED}.gz"
    tabix -f -p bed "${OUT_BED}.gz"
fi

# Check against the BED actually used in the benchmark
if [[ -n "${COMPARE:-}" ]]; then
    if cmp -s "$OUT_BED" "$COMPARE"; then
        echo "MATCH: identical to ${COMPARE}"
    elif diff -q <(sort "$OUT_BED") <(sort "$COMPARE") >/dev/null; then
        echo "SAME REGIONS, DIFFERENT ORDER than ${COMPARE}"
        exit 2
    else
        echo "DIFFERENT regions from ${COMPARE}:"
        diff <(sort "$OUT_BED") <(sort "$COMPARE") | head -20
        exit 1
    fi
fi
