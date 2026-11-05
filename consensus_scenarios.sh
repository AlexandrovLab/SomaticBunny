#!/bin/bash
# Post-EVC consensus calling with six selectable scenarios.
#
# Usage:
#   consensus_scenarios.sh [-s SCENARIOS] [-t TLOD] [-e SOMATIC_EVS] [-p] \
#       SAMPLE MUTECT2_VCF MUSE_VCF STRELKA_SNV STRELKA_INDEL SAGE_VCF \
#       TUMOR_BAM DATABASE_PATH REF_FASTA
#
#   -s  Scenario(s) to run: a number 1-6, a comma list (e.g. 3,5,6) or "all".
#       Default: 6
#   -t  Minimum Mutect2 TLOD for scenario 6 (default 10)
#   -e  Minimum Strelka2 SomaticEVS for scenario 6 (default 13)
#   -p  Scenario 6 drops variants that Mutect2 flagged panel_of_normals.
#       Default: off (PON sites are kept; the flag stays visible in caller_filters)
#
# Scenarios (callers: mt=Mutect2, sa=SAGE, st=Strelka2, ms=MuSE; MuSE is SNV-only,
# so indels have 3 callers):
#   S1  All raw calls from every caller, no filtering (union).
#   S2  Only PASS calls from every caller (union of PASS calls).
#   S3  Called by >= 2 callers, any FILTER status.
#   S4  Called by >= 2 callers, at least one of which is PASS.
#   S5  PASS in >= 2 callers.
#   S6  S5, then (SNVs only) DKFZ bias filter, then TLOD / SomaticEVS thresholds with
#       rescue, then (only with -p) Mutect2 panel-of-normals veto.
#       Rescue: a PASS caller that fails its threshold is demoted; the variant
#       still passes if >= 2 PASS callers remain.
#
# All scenarios: inputs are normalized with bcftools norm (left-align, split
# multiallelics), MNVs are split into SNVs, each caller votes at most once per
# variant (CHROM/POS/REF/ALT), and only primary chromosomes (chr1-22, X, Y, with or without the chr prefix) are kept.
# Outputs, per scenario N:
#   scenarioN/<SAMPLE>_snv.vcf     retained SNVs (valid sites-only VCF)
#   scenarioN/<SAMPLE>_indel.vcf   retained indels
#   scenarioN/<SAMPLE>_snv_annotated.tsv / _indel_annotated.tsv
#       every variant seen by any caller, the scenario decision (PASS or the reason
#       it was dropped), and each caller's FILTER and INFO.
#   special_mutations/             multiallelic / non-SNV records set aside
#   intermediate/                  normalized and split per-caller VCFs

set -euo pipefail

usage() {
    sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//' >&2
    exit 1
}

######################
## ACCEPT ARGUMENTS ##
######################
SCENARIOS=6
TLOD=10
SEVS=13
PON_VETO=0
while getopts ":s:t:e:ph" opt; do
    case $opt in
        s) SCENARIOS=$OPTARG ;;
        t) TLOD=$OPTARG ;;
        e) SEVS=$OPTARG ;;
        p) PON_VETO=1 ;;
        h|*) usage ;;
    esac
done
shift $((OPTIND - 1))
[[ $# -eq 9 ]] || usage

SAMPLE_NAME=$1
MUTECT2_VCF=$2
MUSE_VCF=$3
STRELKA_SNV=$4
STRELKA_INDEL=$5
SAGE_VCF=$6
TUMOR_BAM=$7
DATABASE_PATH=$8
REF_FASTA=$9
PYTHON=${PYTHON:-python}
S=$SAMPLE_NAME

[[ $SCENARIOS == all ]] && SCENARIOS=1,2,3,4,5,6
IFS=',' read -r -a SCENARIO_LIST <<< "$SCENARIOS"
RUN_S6=0
for sc in "${SCENARIO_LIST[@]}"; do
    [[ $sc =~ ^[1-6]$ ]] || { echo "ERROR: invalid scenario '$sc' (use 1-6, a comma list, or all)" >&2; exit 1; }
    [[ $sc == 6 ]] && RUN_S6=1
done
[[ $TLOD =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "ERROR: -t must be a number" >&2; exit 1; }
[[ $SEVS =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "ERROR: -e must be a number" >&2; exit 1; }

# Fail before doing any work if a tool or input is missing.
required_tools=(bcftools awk sort)
required_files=("$MUTECT2_VCF" "$MUSE_VCF" "$STRELKA_SNV" "$STRELKA_INDEL" "$SAGE_VCF" "$REF_FASTA")
if [[ $RUN_S6 -eq 1 ]]; then
    required_tools+=("$PYTHON")
    required_files+=("$TUMOR_BAM" "$DATABASE_PATH/bseq_header" "$DATABASE_PATH/DKFZBiasFilter/scripts/biasFilter.py")
fi
for tool in "${required_tools[@]}"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required tool not found: $tool" >&2; exit 1; }
done
for path in "${required_files[@]}"; do
    [[ -r "$path" ]] || { echo "ERROR: required input not found or unreadable: $path" >&2; exit 1; }
done

WORK_DIR=$PWD
INT_DIR=$WORK_DIR/intermediate
SPECIAL_DIR=$WORK_DIR/special_mutations
mkdir -p "$INT_DIR" "$SPECIAL_DIR"

######################
## Helper functions ##
######################

# Left-align indels and split multiallelic sites; plain-text output.
normalize() {
    bcftools norm -f "$REF_FASTA" -m -any -Ov -o "$2" "$1"
}

# Split a VCF into SNV and indel files. Multiallelic records go to $4.
# Equal-length multi-base substitutions (MNVs) are split into single-base SNVs,
# skipping positions where REF and ALT are the same base.
split_variants() {
    local in=$1 snv_out=$2 indel_out=$3 multi_out=$4
    : > "$snv_out"; : > "$indel_out"; : > "$multi_out"
    awk -F'\t' -v OFS='\t' -v snv="$snv_out" -v indel="$indel_out" -v multi="$multi_out" '
        /^#/                      { print > snv; print > indel; next }
        $5 ~ /,/                  { print > multi; next }
        length($4) != length($5)  { print > indel; next }
        length($4) == 1           { print > snv; next }
        {
            ref = $4; alt = $5; pos = $2
            for (i = 1; i <= length(ref); i++) {
                r = substr(ref, i, 1); a = substr(alt, i, 1)
                if (r == a) continue
                $2 = pos + i - 1; $4 = r; $5 = a
                print > snv
            }
        }
    ' "$in"
}

# Print one caller's records in long format: CHROM POS REF ALT TAG FILTER INFO.
# A FILTER of "." (filters not applied) is recorded as "MISSING" so it is never
# mistaken for "caller absent" or PASS.
emit_long() {
    local tag=$1 file=$2
    awk -F'\t' -v OFS='\t' -v tag="$tag" '
        !/^#/ && $5 != "." && $1 ~ /^(chr)?([0-9]+|X|Y)$/ {
            print $1, $2, $4, $5, tag, ($7 == "." ? "MISSING" : $7), $8
        }
    ' "$file"
}

# Build the per-variant table from long-format records on stdin.
# Columns: 1 CHROM 2 POS 3 ID 4 REF 5 ALT 6 DECISION
#          7-10  FILTER for mt sa st ms ("." = caller did not report it)
#          11-14 INFO   for mt sa st ms
# If one caller has several records for the same variant (e.g. an MNV split plus
# a separate SNV), it still counts once; a PASS record is preferred.
aggregate() {
    awk -F'\t' -v OFS='\t' '
        BEGIN { ntag = split("mt sa st ms", T, " ") }
        {
            k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            if (!(k in seen)) { seen[k] = 1; keys[++nk] = k }
            c = k SUBSEP $5
            if (!(c in F) || (F[c] != "PASS" && $6 == "PASS")) { F[c] = $6; I[c] = $7 }
        }
        END {
            for (n = 1; n <= nk; n++) {
                split(keys[n], f, SUBSEP)
                line = f[1] OFS f[2] OFS "." OFS f[3] OFS f[4] OFS "."
                for (j = 1; j <= ntag; j++) { c = keys[n] SUBSEP T[j]; line = line OFS ((c in F) ? F[c] : ".") }
                for (j = 1; j <= ntag; j++) { c = keys[n] SUBSEP T[j]; line = line OFS ((c in I) ? I[c] : ".") }
                print line
            }
        }
    ' | sort -t$'\t' -k1,1V -k2,2n -k4,4 -k5,5
}

# Apply a scenario rule to the table: column 6 becomes PASS or the drop reason.
select_scenario() {
    local sc=$1 in=$2 out=$3
    awk -F'\t' -v OFS='\t' -v sc="$sc" '
        {
            n = 0; p = 0
            for (c = 7; c <= 10; c++) { if ($c != ".") n++; if ($c == "PASS") p++ }
            d = "PASS"
            if      (sc == 2 && p < 1)          d = "noPASScaller"
            else if (sc == 3 && n < 2)          d = "lt2callers"
            else if (sc == 4 && n < 2)          d = "lt2callers"
            else if (sc == 4 && p < 1)          d = "noPASScaller"
            else if (sc >= 5 && p < 2)          d = "lt2PASScallers"
            $6 = d
            print
        }
    ' "$in" > "$out"
}

# Scenario 6, step 1 (SNVs): run the DKFZ bias filter on the variants that
# passed so far and copy its FILTER verdict into column 6.
apply_bias_filter() {
    local in=$1 out=$2 tmp=$3
    mkdir -p "$tmp" "$tmp/bseq_tmp"   # biasFilter.py writes temp files here and does not create it
    awk -F'\t' -v OFS='\t' '$6 == "PASS" { print $1, $2, ".", $4, $5, ".", "PASS", ".", "." }' "$in" \
        > "$tmp/candidates.vcf"
    local n_before n_after
    n_before=$(wc -l < "$tmp/candidates.vcf")
    echo "  S6 SNVs before bias filter: $n_before"

    if [[ $n_before -eq 0 ]]; then
        echo "  WARNING: no candidate SNVs; skipping bias filter." >&2
        : > "$tmp/bseq.vcf"
    else
        cat "$DATABASE_PATH/bseq_header" "$tmp/candidates.vcf" > "$tmp/candidates.with_header.vcf"
        # With set -e, the script stops here if biasFilter.py fails.
        "$PYTHON" "$DATABASE_PATH/DKFZBiasFilter/scripts/biasFilter.py" \
            "$tmp/candidates.with_header.vcf" \
            "$(readlink -f "$TUMOR_BAM")" \
            "$REF_FASTA" \
            "$tmp/bseq.with_header.vcf" \
            --tempFolder="$tmp/bseq_tmp" \
            --mapq=1 \
            --baseq=1
        awk '!/^#/' "$tmp/bseq.with_header.vcf" > "$tmp/bseq.vcf"
    fi
    n_after=$(awk -F'\t' '$7 == "PASS"' "$tmp/bseq.vcf" | wc -l)
    echo "  S6 SNVs PASS after bias filter: $n_after"

    # Candidates the bias filter did not return at all are marked bseqMissing.
    awk -F'\t' -v OFS='\t' '
        FILENAME == ARGV[1] { b[$1 SUBSEP $2 SUBSEP $4 SUBSEP $5] = $7; next }
        $6 == "PASS" {
            k = $1 SUBSEP $2 SUBSEP $4 SUBSEP $5
            $6 = (k in b) ? b[k] : "bseqMissing"
        }
        { print }
    ' "$tmp/bseq.vcf" "$in" > "$out"
    rm -rf "$tmp"
}

# Scenario 6, step 2: TLOD / SomaticEVS thresholds with rescue, then (with -p) PON veto.
# Only INFO from PASS callers is checked. A PASS caller below a threshold is
# demoted (its FILTER column becomes lowTLOD / lowSomaticEVS). If >= 2 PASS
# callers remain the variant keeps its status; otherwise the flags are added
# to the variant. With -p, Mutect2 panel_of_normals vetoes the variant.
apply_thresholds() {
    local in=$1 out=$2
    awk -F'\t' -v OFS='\t' -v tlod="$TLOD" -v sevs="$SEVS" -v pon="$PON_VETO" '
        function addflag(f) { if ($6 !~ ("(^|;)" f "(;|$)")) $6 = ($6 == "PASS" ? f : $6 ";" f) }
        $6 == "lt2PASScallers" { print; next }
        {
            npass = 0; nlow = 0; delete low
            for (j = 0; j < 4; j++) {
                fc = 7 + j; ic = 11 + j
                if ($fc != "PASS") continue
                npass++
                reason = ""
                n = split($ic, kv, ";")
                for (i = 1; i <= n; i++) {
                    split(kv[i], p, "=")
                    if (p[1] == "TLOD"       && p[2] + 0 < tlod + 0) reason = (reason == "" ? "" : reason ";") "lowTLOD"
                    if (p[1] == "SomaticEVS" && p[2] + 0 < sevs + 0) reason = (reason == "" ? "" : reason ";") "lowSomaticEVS"
                }
                if (reason != "") { low[j] = reason; nlow++ }
            }
            if (nlow > 0) {
                if (npass - nlow >= 2) {
                    for (j in low) $(7 + j) = low[j]          # rescue: demote that caller
                } else {
                    for (j in low) { m = split(low[j], r, ";"); for (i = 1; i <= m; i++) addflag(r[i]) }
                }
            }
            if (pon == 1 && $7 ~ /(^|;)panel_of_normals(;|$)/) addflag("panel_of_normals")
            print
        }
    ' "$in" > "$out"
}

# Write the annotated TSV and the sites-only VCF of retained variants.
# VCF FILTER uses the same three values in every scenario, based on how many
# callers called the variant PASS (scenario 6: after threshold demotion):
#   PASS           >= 2 PASS callers
#   OnePassCaller  exactly 1 PASS caller
#   NoPassCaller   no PASS caller
# So S5/S6 files contain only PASS, while S1-S4 files show at a glance which
# records would also meet the PASS consensus. INFO holds the summary fields
# followed by the full INFO of every caller that reported the variant, with each
# key prefixed by the caller tag (mt_TLOD, sa_DP, st_SomaticEVS, ms_SOMATIC ...)
# so keys from different callers never clash. The prefixed keys are declared in
# the header from each caller's own ##INFO lines; a key a caller used without
# declaring it is declared generically (Flag, or String with Number=.).
# Args: scenario table tsv vcf mt_file sa_file st_file ms_file
# (pass /dev/null for a caller that has no file, e.g. MuSE for indels).
write_outputs() {
    local sc=$1 table=$2 tsv=$3 vcf=$4
    shift 4
    {
        printf '#CHROM\tPOS\tID\tREF\tALT\tSCENARIO%s_DECISION\tMUTECT2_FILTER\tSAGE_FILTER\tSTRELKA_FILTER\tMUSE_FILTER\tMUTECT2_INFO\tSAGE_INFO\tSTRELKA_INFO\tMUSE_INFO\n' "$sc"
        cat "$table"
    } > "$tsv"

    awk -F'\t' -v OFS='\t' -v sc="$sc" -v pon="$PON_VETO" -v ref="$(readlink -f "$REF_FASTA")" -v contigs="$INT_DIR/${S}_mutect.vcf" '
        BEGIN {
            split("mt sa st ms", T, " ")
            split("Mutect2 SAGE Strelka2 MuSE", NAME, " ")
            ncf = 4   # the first four files are the caller VCFs (headers only are read)
        }
        # Caller headers: turn ##INFO=<ID=X,... into ##INFO=<ID=tag_X,... "[Caller] ...".
        FILENAME != ARGV[ARGC - 1] {
            if ($0 !~ /^#/) next
            if ($0 ~ /^##INFO=<ID=/) {
                for (j = 1; j <= ncf; j++) if (FILENAME == ARGV[j]) break
                line = $0
                id = line; sub(/^##INFO=<ID=/, "", id); sub(/,.*/, "", id)
                key = T[j] "_" id
                if (key in decl) next
                sub(/^##INFO=<ID=[^,]*/, "##INFO=<ID=" key, line)
                sub(/Description="/, "Description=\"[" NAME[j] "] ", line)
                decl[key] = line; order[++nd] = key
            }
            next
        }
        # Table rows: keep retained variants and build their INFO.
        $6 == "PASS" {
            cl = ""; pc = ""; cf = ""; extra = ""; n = 0; p = 0
            for (j = 1; j <= 4; j++) {
                f = $(6 + j)
                if (f == ".") continue
                n++; cl = cl (cl == "" ? "" : ",") T[j]
                g = f; gsub(";", "|", g); cf = cf (cf == "" ? "" : ",") T[j] ":" g
                if (f == "PASS") { p++; pc = pc (pc == "" ? "" : ",") T[j] }
                info = $(10 + j)
                if (info == ".") continue
                m = split(info, kv, ";")
                for (i = 1; i <= m; i++) {
                    if (kv[i] == "" || kv[i] == ".") continue
                    eq = index(kv[i], "=")
                    k = (eq ? substr(kv[i], 1, eq - 1) : kv[i])
                    key = T[j] "_" k
                    if (!(key in decl) && !(key in undecl)) {
                        undecl[key] = sprintf("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"[%s] %s (not declared in caller header)\">", \
                            key, (eq ? "." : "0"), (eq ? "String" : "Flag"), NAME[j], k)
                        uorder[++nu] = key
                    }
                    extra = extra ";" T[j] "_" kv[i]
                }
            }
            out = "callers=" cl ";ncallers=" n ";npass=" p
            if (p > 0) out = out ";pass_callers=" pc
            out = out ";caller_filters=" cf extra
            filt = (p >= 2 ? "PASS" : (p == 1 ? "OnePassCaller" : "NoPassCaller"))
            rec[++nr] = $1 OFS $2 OFS "." OFS $4 OFS $5 OFS "." OFS filt OFS out
        }
        END {
            print "##fileformat=VCFv4.2"
            print "##source=consensus_scenarios.sh scenario " sc
            print "##reference=file://" ref
            while ((getline line < contigs) > 0) { if (line ~ /^##contig=/) print line; else if (line !~ /^##/) break }
            print "##FILTER=<ID=PASS,Description=\"Called PASS by >= 2 callers" (sc == 6 ? "; also passed scenario 6 bias filter and TLOD/SomaticEVS checks" (pon == 1 ? " and the Mutect2 panel-of-normals veto" : "") : "") "\">"
            print "##FILTER=<ID=OnePassCaller,Description=\"Called PASS by exactly 1 caller (see caller_filters)\">"
            print "##FILTER=<ID=NoPassCaller,Description=\"Not called PASS by any caller (see caller_filters)\">"
            print "##INFO=<ID=callers,Number=.,Type=String,Description=\"Callers reporting the variant, any FILTER (mt=Mutect2, sa=SAGE, st=Strelka2, ms=MuSE)\">"
            print "##INFO=<ID=ncallers,Number=1,Type=Integer,Description=\"Number of distinct callers reporting the variant\">"
            print "##INFO=<ID=pass_callers,Number=.,Type=String,Description=\"Callers with FILTER=PASS (scenario 6: after threshold demotion)\">"
            print "##INFO=<ID=npass,Number=1,Type=Integer,Description=\"Number of distinct PASS callers\">"
            print "##INFO=<ID=caller_filters,Number=.,Type=String,Description=\"Per-caller FILTER as caller:filter; multiple filters joined with |\">"
            for (i = 1; i <= nd; i++) print decl[order[i]]
            for (i = 1; i <= nu; i++) print undecl[uorder[i]]
            print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            for (i = 1; i <= nr; i++) print rec[i]
        }
    ' "$@" "$table" > "$vcf"
}

#######################
## Normalize inputs ##
#######################
echo "Normalizing input VCFs..."
normalize "$MUTECT2_VCF"   "$INT_DIR/${S}_mutect.vcf"
normalize "$MUSE_VCF"      "$INT_DIR/${S}_muse_snv.vcf"
normalize "$STRELKA_SNV"   "$INT_DIR/${S}_strelka_snv.vcf"
normalize "$STRELKA_INDEL" "$INT_DIR/${S}_strelka_indel.vcf"
normalize "$SAGE_VCF"      "$INT_DIR/${S}_sage.vcf"

#######################
## Fix Mutect2 files ##
#######################
# Put the normal sample in column 10 and the tumor in column 11, then split SNVs/indels.
echo "Fixing Mutect2 files..."
f=$INT_DIR/${S}_mutect.vcf
tumor_name=$(awk -F= '/^##tumor_sample=/ { print $2; exit }' "$f")
first_sample=$(awk -F'\t' '/^#CHROM/ { print $10; exit }' "$f")
if [[ -n "$tumor_name" && "$first_sample" == "$tumor_name" ]]; then
    echo "  Wrong order. Swapping columns..."
    swap=1
else
    echo "  Correct order (or tumor name not found in header). Not swapping."
    swap=0
fi
awk -F'\t' -v OFS='\t' -v swap="$swap" '
    /^##/ { print; next }
    swap  { t = $10; $10 = $11; $11 = t }
    { print }
' "$f" > "$INT_DIR/${S}_mutect.swapped.vcf"
split_variants "$INT_DIR/${S}_mutect.swapped.vcf" \
    "$INT_DIR/${S}_mutect_snv.vcf" "$INT_DIR/${S}_mutect_indel.vcf" "$SPECIAL_DIR/${S}_mutect.multiMutation.txt"

####################
## Fix MuSE files ##
####################
# MuSE writes tumor first: swap to normal/tumor and set ID to ".".
# Anything that is not a single-base SNV is set aside, never truncated.
echo "Fixing MuSE files..."
f=$INT_DIR/${S}_muse_snv.vcf
: > "$SPECIAL_DIR/${S}_muse_snvMulti.txt"
awk -F'\t' -v OFS='\t' -v multi="$SPECIAL_DIR/${S}_muse_snvMulti.txt" '
    /^##/     { print; next }
    /^#CHROM/ { t = $10; $10 = $11; $11 = t; print; next }
    {
        if (length($4) != 1 || length($5) != 1 || $4 == $5) { print > multi; next }
        $3 = "."
        t = $10; $10 = $11; $11 = t
        print
    }
' "$f" > "$f.tmp"
mv "$f.tmp" "$f"

####################
## Fix SAGE files ##
####################
echo "Fixing SAGE files..."
split_variants "$INT_DIR/${S}_sage.vcf" \
    "$INT_DIR/${S}_sage_snv.vcf" "$INT_DIR/${S}_sage_indel.vcf" "$SPECIAL_DIR/${S}_sage.multiMutation.txt"

######################################
## Per-caller summary (informational) ##
######################################
echo "=== Per-caller records / PASS ==="
for f in "$INT_DIR/${S}_mutect_snv.vcf" "$INT_DIR/${S}_muse_snv.vcf" "$INT_DIR/${S}_strelka_snv.vcf" "$INT_DIR/${S}_sage_snv.vcf" \
         "$INT_DIR/${S}_mutect_indel.vcf" "$INT_DIR/${S}_strelka_indel.vcf" "$INT_DIR/${S}_sage_indel.vcf"; do
    awk -F'\t' -v name="$(basename "$f")" '!/^#/ { n++; if ($7 == "PASS") p++ }
        END { printf "  %s: %d records, %d PASS\n", name, n, p }' "$f"
done

#############################
## Build per-variant tables ##
#############################
echo "Building per-variant tables..."
{
    emit_long mt "$INT_DIR/${S}_mutect_snv.vcf"
    emit_long sa "$INT_DIR/${S}_sage_snv.vcf"
    emit_long st "$INT_DIR/${S}_strelka_snv.vcf"
    emit_long ms "$INT_DIR/${S}_muse_snv.vcf"
} | aggregate > "$INT_DIR/${S}_snv_table.tsv"

{
    emit_long mt "$INT_DIR/${S}_mutect_indel.vcf"
    emit_long sa "$INT_DIR/${S}_sage_indel.vcf"
    emit_long st "$INT_DIR/${S}_strelka_indel.vcf"
} | aggregate > "$INT_DIR/${S}_indel_table.tsv"

echo "  Distinct SNVs: $(wc -l < "$INT_DIR/${S}_snv_table.tsv"), distinct indels: $(wc -l < "$INT_DIR/${S}_indel_table.tsv")"

###################
## Run scenarios ##
###################
for sc in "${SCENARIO_LIST[@]}"; do
    echo "Running scenario $sc..."
    OUT=$WORK_DIR/scenario$sc
    mkdir -p "$OUT"

    for vt in snv indel; do
        select_scenario "$sc" "$INT_DIR/${S}_${vt}_table.tsv" "$OUT/.${vt}.step1.tsv"
        if [[ $sc -eq 6 ]]; then
            if [[ $vt == snv ]]; then
                apply_bias_filter "$OUT/.${vt}.step1.tsv" "$OUT/.${vt}.step2.tsv" "$OUT/.bseq_work"
            else
                cp "$OUT/.${vt}.step1.tsv" "$OUT/.${vt}.step2.tsv"   # bias filter is SNV-only
            fi
            apply_thresholds "$OUT/.${vt}.step2.tsv" "$OUT/.${vt}.final.tsv"
        else
            cp "$OUT/.${vt}.step1.tsv" "$OUT/.${vt}.final.tsv"
        fi
        if [[ $vt == snv ]]; then
            caller_files=("$INT_DIR/${S}_mutect_snv.vcf" "$INT_DIR/${S}_sage_snv.vcf" "$INT_DIR/${S}_strelka_snv.vcf" "$INT_DIR/${S}_muse_snv.vcf")
        else
            caller_files=("$INT_DIR/${S}_mutect_indel.vcf" "$INT_DIR/${S}_sage_indel.vcf" "$INT_DIR/${S}_strelka_indel.vcf" /dev/null)
        fi
        write_outputs "$sc" "$OUT/.${vt}.final.tsv" "$OUT/${S}_${vt}_annotated.tsv" "$OUT/${S}_${vt}.vcf" "${caller_files[@]}"
        rm -f "$OUT"/.${vt}.*.tsv
        echo "  ${vt}s retained: $(awk '!/^#/' "$OUT/${S}_${vt}.vcf" | wc -l)"
    done
done

echo "Processing complete for ${S}"
