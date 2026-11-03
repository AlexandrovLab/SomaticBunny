#!/bin/bash
# Run once on TSCC from the repo root: whole-contig ApplyBQSR interval lists for the packaged references.
# Writes <interval_dir>/bqsr/<k>-contigs.interval_list, k = 1..20, next to the existing <k>-scattered.interval_list files.
set -euo pipefail
DB=${1:?usage: make_packaged_intervals.sh /path/to/SomaticBunny_database}
mk() { echo "== $1"; python3 bin/make_contig_groups.py "$2" 20 "$DB/Databases/$1/$3/bqsr"; }
mk GRCh38 "$DB/GRCh38_ref/GRCh38.d1.vd1.dict" GRCh38_interval_list_20
mk GRCh37 "$DB/GRCh37_ref/GRCh37.dict"        GRCh37_interval_list_20
mk mm39   "$DB/mm39_ref/mm39.dict"            mm39_interval_list_20
mk RN7    "$DB/RN7_ref/rn7.dict"              RN7_interval_list_20
