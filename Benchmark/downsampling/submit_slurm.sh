#!/bin/bash
#SBATCH --job-name=downsample_wgs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 80:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/downsample_%j.out
#SBATCH --error=logs/downsample_%j.err
# Adjust for your cluster, e.g.:
# #SBATCH -p <partition>
# #SBATCH -A <account>

set -euo pipefail

SAMPLE_LIST="samples.txt"
INPUT_DIR="/path/to/original_fastq"
OUTPUT_DIR="/path/to/fastq_downsampled_60x_or_30x"

# export SEQTK=/path/to/seqtk   

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"

bash "${SCRIPT_DIR}/downsample.sh" \
    -s "$SAMPLE_LIST" \
    -i "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    --jobs 4
