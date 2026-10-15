#!/usr/bin/env bash
#######################################################################
# SomaticBunny benchmark: SomaticBunny on WES
#
# NOTE: verify the exome-mode options against the actual benchmark
# session (--type exome and the WES_BED paths in paths.env are
# assumed from the WGS command; adjust if the WES run differed).
#######################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOMATICBUNNY_DIR="${SOMATICBUNNY_DIR:?set SOMATICBUNNY_DIR to the SomaticBunny checkout (directory containing main.nf)}"
PATHS_ENV="${PATHS_ENV:-${SCRIPT_DIR}/paths.env}"

# ---- Nextflow environment (benchmark: conda env "env_nf") ----
if [[ -z "${CONDA_SH:-}" ]]; then
    CONDA_SH="$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh"
fi
# shellcheck disable=SC1090
source "$CONDA_SH" || { echo "ERROR: cannot source conda: $CONDA_SH (set CONDA_SH=...)" >&2; exit 1; }
conda activate "${NF_CONDA_ENV:-env_nf}" || { echo "ERROR: cannot activate conda env ${NF_CONDA_ENV:-env_nf}" >&2; exit 1; }

set -euo pipefail

export NXF_OPTS='-Xms1g -Xmx4g'

# ---- reference paths (exports TMPDIR, as in the original session) ----
[[ -f "$PATHS_ENV" ]] || { echo "ERROR: paths file not found: $PATHS_ENV" >&2; exit 1; }
set -a
# shellcheck disable=SC1090
source "$PATHS_ENV"
set +a
mkdir -p "$TMPDIR" 2>/dev/null || { echo "ERROR: cannot create TMPDIR=${TMPDIR} (edit paths.env)" >&2; exit 1; }

missing=0
for f in "${SOMATICBUNNY_DIR}/main.nf" \
         "$REF" "$REF_FAI" "$REF_DICT" \
         "$WES_BED" "$WES_BED_TBI"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; missing=1; }
done
[[ -d "$DATABASE_PATH" ]] || { echo "ERROR: database path not found: $DATABASE_PATH" >&2; missing=1; }
(( missing == 0 )) || exit 1

command -v nextflow >/dev/null 2>&1 || { echo "ERROR: nextflow not on PATH" >&2; exit 1; }
nextflow -version 2>&1 | grep -i "version" | head -n1 > nextflow_version.txt || true

nextflow run "${SOMATICBUNNY_DIR}/main.nf" -profile tscc_slurm \
    --type exome --first_step mapping --genome GRCh38 \
    --publish recal_bam \
    --database_path "$DATABASE_PATH" \
    --ref "$REF" \
    --ref_fai "$REF_FAI" \
    --ref_dict "$REF_DICT" \
    --bed "$WES_BED" \
    --bed_tbi "$WES_BED_TBI" \
    --tool ascat,manta,delly,cnvkit \
    -with-report \
    "$@"
