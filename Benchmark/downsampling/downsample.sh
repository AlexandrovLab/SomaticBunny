#!/usr/bin/env bash
#######################################################################
# SomaticBunny benchmark: downsample paired-end tumor/normal WGS FASTQs
#
# Randomly subsamples read pairs with `seqtk sample` to reach a target
# coverage. The same seed is used for R1 and R2 so mates stay paired.
#
# Defaults reproduce the published benchmark:
#   tumor  90x -> 60x  (fraction .6666)
#   normal 40x -> 30x  (fraction .7500)
#   seed   1234
#
# Usage: see ./downsample.sh --help
#######################################################################
set -uo pipefail
ORIG_CMD="$0 $*"

# ---------------- defaults (values used in the benchmark) ----------------
SEED=1234
TUMOR_CURRENT=90
TUMOR_TARGET=60
NORMAL_CURRENT=40
NORMAL_TARGET=30
MAX_JOBS=4
WHICH="both"          # both | tumor | normal
SEQTK="${SEQTK:-seqtk}" 

SAMPLE_LIST=""
INPUT_DIR=""
OUTPUT_DIR=""

usage() {
    cat <<USAGE
Usage: $(basename "$0") -s SAMPLE_LIST -i INPUT_DIR -o OUTPUT_DIR [options]

Required:
  -s, --samples FILE        One sample ID per line (lines starting with # are ignored)
  -i, --input-dir DIR       Directory with {sample}_{tumor,normal}_{1,2}.fastq.gz
  -o, --output-dir DIR      Directory for downsampled FASTQs, logs and run info

Options:
  --seed INT                Random seed for seqtk (default: ${SEED})
  --tumor-current NUM       Tumor input coverage (default: ${TUMOR_CURRENT})
  --tumor-target NUM        Tumor target coverage (default: ${TUMOR_TARGET})
  --normal-current NUM      Normal input coverage (default: ${NORMAL_CURRENT})
  --normal-target NUM       Normal target coverage (default: ${NORMAL_TARGET})
  --which both|tumor|normal Which samples to downsample (default: ${WHICH})
  -j, --jobs INT            Samples processed in parallel (default: ${MAX_JOBS})
  -h, --help                Show this help

Environment:
  SEQTK                     Path to the seqtk binary (default: seqtk on PATH)
USAGE
}

# ---------------- argument parsing ----------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samples)       SAMPLE_LIST="$2"; shift 2 ;;
        -i|--input-dir)     INPUT_DIR="$2"; shift 2 ;;
        -o|--output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --seed)             SEED="$2"; shift 2 ;;
        --tumor-current)    TUMOR_CURRENT="$2"; shift 2 ;;
        --tumor-target)     TUMOR_TARGET="$2"; shift 2 ;;
        --normal-current)   NORMAL_CURRENT="$2"; shift 2 ;;
        --normal-target)    NORMAL_TARGET="$2"; shift 2 ;;
        --which)            WHICH="$2"; shift 2 ;;
        -j|--jobs)          MAX_JOBS="$2"; shift 2 ;;
        -h|--help)          usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# ---------------- checks ----------------
[[ -z "$SAMPLE_LIST" || -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]] && { usage >&2; exit 1; }
[[ -f "$SAMPLE_LIST" ]] || { echo "ERROR: sample list not found: $SAMPLE_LIST" >&2; exit 1; }
[[ -d "$INPUT_DIR" ]]   || { echo "ERROR: input dir not found: $INPUT_DIR" >&2; exit 1; }
[[ "$WHICH" =~ ^(both|tumor|normal)$ ]] || { echo "ERROR: --which must be both, tumor or normal" >&2; exit 1; }
command -v "$SEQTK" >/dev/null 2>&1 || { echo "ERROR: seqtk not found (set SEQTK=/path/to/seqtk)" >&2; exit 1; }
command -v bc >/dev/null 2>&1 || { echo "ERROR: bc not found" >&2; exit 1; }

# Kept identical to the original benchmark run so outputs are reproducible.
TUMOR_FRACTION=$(echo "scale=4; $TUMOR_TARGET / $TUMOR_CURRENT" | bc)
NORMAL_FRACTION=$(echo "scale=4; $NORMAL_TARGET / $NORMAL_CURRENT" | bc)

SEQTK_VERSION=$("$SEQTK" 2>&1 | awk -F': ' '/^Version/{print $2}')
ts() { date '+%Y-%m-%d %H:%M:%S'; }

mkdir -p "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/manifest"

# ---------------- record provenance ----------------
{
    echo "date             : $(ts)"
    echo "host             : $(hostname)"
    echo "command          : ${ORIG_CMD}"
    echo "seqtk            : ${SEQTK} (${SEQTK_VERSION:-unknown version})"
    echo "gzip             : $(gzip --version 2>&1 | head -n1)"
    echo "sample_list      : ${SAMPLE_LIST}"
    echo "sample_list_md5  : $(md5sum "$SAMPLE_LIST" | cut -d' ' -f1)"
    echo "input_dir        : ${INPUT_DIR}"
    echo "output_dir       : ${OUTPUT_DIR}"
    echo "seed             : ${SEED}"
    echo "which            : ${WHICH}"
    echo "tumor            : ${TUMOR_CURRENT}x -> ${TUMOR_TARGET}x (fraction ${TUMOR_FRACTION})"
    echo "normal           : ${NORMAL_CURRENT}x -> ${NORMAL_TARGET}x (fraction ${NORMAL_FRACTION})"
} > "${OUTPUT_DIR}/run_info.txt"

echo "=========================================="
echo "SomaticBunny benchmark: WGS downsampling"
echo "=========================================="
cat "${OUTPUT_DIR}/run_info.txt"
echo "parallel jobs    : ${MAX_JOBS}"
echo "=========================================="

# ---------------- downsample one tumor or normal pair ----------------
downsample_pair() {
    local sample=$1 status=$2 fraction=$3
    local read in tmp nlines
    local -a n=()

    for read in 1 2; do
        in="${INPUT_DIR}/${sample}_${status}_${read}.fastq.gz"
        [[ -f "$in" ]] || { echo "  ERROR: missing input $in"; return 1; }
    done

    echo "  [$(ts)] Downsampling ${status} (fraction ${fraction}, seed ${SEED})"
    for read in 1 2; do
        in="${INPUT_DIR}/${sample}_${status}_${read}.fastq.gz"
        tmp="${OUTPUT_DIR}/${sample}_${status}_${read}.fastq.gz.tmp"
        if ! gzip -cd "$in" \
            | "$SEQTK" sample -s "$SEED" - "$fraction" \
            | awk -v cf="${tmp}.lines" '{ print } END { print NR > cf }' \
            | gzip > "$tmp"; then
            echo "  ERROR: failed on $in (corrupt/truncated input, or write error)"
            rm -f "${OUTPUT_DIR}/${sample}_${status}"_[12].fastq.gz.tmp{,.lines}
            return 1
        fi
        nlines=$(cat "${tmp}.lines")
        if (( nlines == 0 || nlines % 4 != 0 )); then
            echo "  ERROR: ${tmp%.tmp} has ${nlines} lines (empty or not FASTQ)"
            rm -f "${OUTPUT_DIR}/${sample}_${status}"_[12].fastq.gz.tmp{,.lines}
            return 1
        fi
        n[$read]=$(( nlines / 4 ))
    done

    if [[ "${n[1]}" != "${n[2]}" ]]; then
        echo "  ERROR: R1 (${n[1]}) and R2 (${n[2]}) read counts differ for ${sample} ${status}"
        rm -f "${OUTPUT_DIR}/${sample}_${status}"_[12].fastq.gz.tmp{,.lines}
        return 1
    fi

    # Publish both files only after both passed
    for read in 1 2; do
        mv -f "${OUTPUT_DIR}/${sample}_${status}_${read}.fastq.gz.tmp" \
              "${OUTPUT_DIR}/${sample}_${status}_${read}.fastq.gz"
        rm -f "${OUTPUT_DIR}/${sample}_${status}_${read}.fastq.gz.tmp.lines"
    done
    echo "  ${n[1]} read pairs"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$status" "$fraction" "$SEED" "${n[1]}" "${n[2]}" \
        > "${OUTPUT_DIR}/manifest/${sample}_${status}.tsv"
    echo "  [$(ts)] ${status} done"
}

# ---------------- process one sample (tumor + normal) ----------------
process_sample() {
    local sample=$1 rc=0
    local log="${OUTPUT_DIR}/logs/${sample}.log"
    {
        echo "[$(ts)] Processing ${sample}"
        if [[ "$WHICH" != "normal" ]]; then
            downsample_pair "$sample" tumor "$TUMOR_FRACTION" || rc=1
        fi
        if [[ "$WHICH" != "tumor" ]]; then
            downsample_pair "$sample" normal "$NORMAL_FRACTION" || rc=1
        fi
        echo "[$(ts)] Finished ${sample} (exit ${rc})"
    } > "$log" 2>&1

    if [[ $rc -eq 0 ]]; then
        echo "[$(ts)] ${sample}: OK"
    else
        echo "[$(ts)] ${sample}: FAILED (see ${log})"
    fi
    return $rc
}

# ---------------- main loop ----------------
pids=()
n_samples=0
while IFS= read -r sample || [[ -n "$sample" ]]; do
    sample="${sample%$'\r'}"                     
    sample="$(echo "$sample" | xargs)"         
    [[ -z "$sample" || "$sample" == \#* ]] && continue

    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 5
    done

    process_sample "$sample" < /dev/null &
    pids+=("$!")
    n_samples=$((n_samples + 1))
done < "$SAMPLE_LIST"

echo ""
echo "Waiting for all jobs to complete..."
failed=0
for pid in ${pids[@]+"${pids[@]}"}; do
    wait "$pid" || failed=$((failed + 1))
done

{
    printf 'sample\tstatus\tfraction\tseed\treads_R1\treads_R2\n'
    cat "${OUTPUT_DIR}"/manifest/*.tsv 2>/dev/null | sort
} > "${OUTPUT_DIR}/manifest.tsv"

echo ""
echo "=========================================="
echo "Samples processed : ${n_samples}"
echo "Samples failed    : ${failed}"
echo "FASTQs written    : $(ls -1 "${OUTPUT_DIR}"/*.fastq.gz 2>/dev/null | wc -l)"
echo "Manifest          : ${OUTPUT_DIR}/manifest.tsv"
echo "Run info          : ${OUTPUT_DIR}/run_info.txt"
echo "Logs              : ${OUTPUT_DIR}/logs"
echo "Done: $(ts)"
echo "=========================================="

[[ $failed -eq 0 ]]
