#!/bin/bash
# Post-EVC consensus calling with six selectable scenarios.
#
# Usage:
#   consensus_scenarios.sh [-s SCENARIOS] [-t TLOD] [-e SOMATIC_EVS] [-p] \
#       [-c CONTIGS] [-b CALLABLE_BED] \
#       SAMPLE MUTECT2_VCF MUSE_VCF STRELKA_SNV STRELKA_INDEL SAGE_VCF \
#       TUMOR_BAM DATABASE_PATH REF_FASTA
#
#   -s  Scenario(s) to run: a number 1-6, a comma list (e.g. 3,5,6) or "all".
#       Default: 6
#   -t  Minimum Mutect2 TLOD for scenario 6 (default 10)
#   -e  Minimum Strelka2 SomaticEVS for scenario 6 (default 13)
#   -p  Scenario 6 drops variants that Mutect2 flagged panel_of_normals.
#       Default: off (PON sites are kept; the flag stays visible in caller_filters)
#   -c  Contigs allowed into the consensus: a comma list (chr1,chr2,...) or a file
#       with one contig name per line. Every name must exist in REF_FASTA.fai.
#   -b  Callable-region BED (0-based, half-open). Only variants whose POS lies in
#       a region enter the consensus. Every BED contig must exist in REF_FASTA.fai.
#       With -c as well, only BED regions on the -c contigs are used.
#       Default (neither -c nor -b): the human primary chromosomes chr1-22, chrX,
#       chrY (or 1-22, X, Y), which must all be present in REF_FASTA.fai; any
#       other reference must be given -c or -b explicitly.
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
# multiallelics) and only variants inside the allowed contigs / callable regions
# (-c / -b) enter the consensus; everything else (viral, EBV, decoy, alt, chrM ...)
# is listed in special_mutations/<SAMPLE>_outside_callable_regions.tsv.
# MNVs are decomposed into SNVs so that SNV-only callers can support them; each
# derived SNV is flagged MNV_DERIVED, keeps the caller's FILTER, and drops the
# caller's event-level INFO and QUAL (TLOD, SomaticEVS ... were computed for the
# whole MNV), so it is score-neutral in scenario 6. The original MNV records are
# kept in special_mutations/<SAMPLE>_<caller>.mnv.vcf.
# Each caller votes at most once per variant (CHROM/POS/REF/ALT). If one caller
# has several records for the same variant, the kept record is chosen by a fixed
# order: PASS first, then a directly called record over an MNV-derived one, then
# the higher caller score (Mutect2 TLOD, Strelka2 SomaticEVS, SAGE/MuSE QUAL),
# then the lexically smallest FILTER/INFO/QUAL text, so the result never depends
# on input order. The number of records per caller is kept in caller_records.
# Outputs, per scenario N:
#   scenarioN/<SAMPLE>_snv.vcf     retained SNVs (valid sites-only VCF)
#   scenarioN/<SAMPLE>_indel.vcf   retained indels
#   scenarioN/<SAMPLE>_snv_annotated.tsv / _indel_annotated.tsv
#       every variant seen by any caller, the scenario decision (PASS or the reason
#       it was dropped), and each caller's FILTER and INFO.
#   special_mutations/             multiallelic / non-SNV records set aside
#   intermediate/                  normalized and split per-caller VCFs

set -euo pipefail
export LC_ALL=C   # byte-order string comparison and sorting, so results never depend on the locale

usage() {
    sed -n '2,59p' "$0" | sed 's/^# \{0,1\}//' >&2
    exit 1
}

######################
## ACCEPT ARGUMENTS ##
######################
SCENARIOS=6
TLOD=10
SEVS=13
PON_VETO=0
CONTIGS_ARG=""
CALLABLE_BED=""
while getopts ":s:t:e:c:b:ph" opt; do
    case $opt in
        s) SCENARIOS=$OPTARG ;;
        t) TLOD=$OPTARG ;;
        e) SEVS=$OPTARG ;;
        p) PON_VETO=1 ;;
        c) CONTIGS_ARG=$OPTARG ;;
        b) CALLABLE_BED=$OPTARG ;;
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

##############################################
## Allowed contigs / callable regions (-c/-b) ##
##############################################
# Everything is checked against the reference index; any mismatch is fatal.
FAI=${REF_FASTA}.fai
[[ -r "$FAI" ]] || { echo "ERROR: reference index not found: $FAI (needed to validate the allowed contigs)" >&2; exit 1; }
ALLOWED_CONTIGS=$INT_DIR/allowed_contigs.txt
REGIONS=$INT_DIR/callable_regions.bed

if [[ -n "$CONTIGS_ARG" ]]; then
    if [[ -f "$CONTIGS_ARG" ]]; then
        awk '{ sub(/#.*/, ""); gsub(/[ \t\r]/, "") } NF' "$CONTIGS_ARG" > "$ALLOWED_CONTIGS"
        REGION_SOURCE="contig list $(readlink -f "$CONTIGS_ARG")"
    else
        tr ',' '\n' <<< "$CONTIGS_ARG" | awk 'NF' > "$ALLOWED_CONTIGS"
        REGION_SOURCE="contigs $CONTIGS_ARG"
    fi
    [[ -s "$ALLOWED_CONTIGS" ]] || { echo "ERROR: -c gave no contig names" >&2; exit 1; }
elif [[ -z "$CALLABLE_BED" ]]; then
    # Default: the human primary chromosomes, using the exact names in the reference.
    for prefix in chr ""; do
        printf "${prefix}%s\n" $(seq 1 22) X Y > "$ALLOWED_CONTIGS"
        if [[ $(awk -F'\t' 'FILENAME == ARGV[1] { want[$1] = 1; next } ($1 in want) { n++ } END { print n + 0 }' \
                "$ALLOWED_CONTIGS" "$FAI") -eq 24 ]]; then
            REGION_SOURCE="human primary chromosomes ${prefix}1-22, ${prefix}X, ${prefix}Y (default)"
            break
        fi
        : > "$ALLOWED_CONTIGS"
    done
    [[ -s "$ALLOWED_CONTIGS" ]] || {
        echo "ERROR: $FAI does not contain the human primary chromosomes (chr1-22, chrX, chrY or 1-22, X, Y)." >&2
        echo "       Give the contigs to analyse with -c, or a callable-region BED with -b." >&2
        exit 1; }
fi

if [[ -s "$ALLOWED_CONTIGS" ]]; then
    missing=$(awk -F'\t' 'FILENAME == ARGV[1] { have[$1] = 1; next } !($1 in have)' "$FAI" "$ALLOWED_CONTIGS" | sort -u | paste -sd, -)
    [[ -z "$missing" ]] || { echo "ERROR: allowed contig(s) not in $FAI: $missing" >&2; exit 1; }
fi

if [[ -n "$CALLABLE_BED" ]]; then
    [[ -r "$CALLABLE_BED" ]] || { echo "ERROR: callable BED not found or unreadable: $CALLABLE_BED" >&2; exit 1; }
    # Validate: 3+ columns, integer coordinates, start < end, end <= contig length, contig in the reference.
    awk -F'\t' -v OFS='\t' -v keep="$ALLOWED_CONTIGS" '
        BEGIN { while ((getline line < keep) > 0) { use[line] = 1; nuse++ } }
        FILENAME == ARGV[1] { len[$1] = $2; next }
        /^(#|track|browser)/ || NF == 0 { next }
        {
            if (NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 + 0 >= $3 + 0) {
                printf "ERROR: malformed BED line %d: %s\n", FNR, $0 > "/dev/stderr"; bad = 1; exit 1 }
            if (!($1 in len)) { printf "ERROR: BED contig not in the reference index: %s\n", $1 > "/dev/stderr"; bad = 1; exit 1 }
            if ($3 + 0 > len[$1] + 0) { printf "ERROR: BED region beyond the end of %s (%s > %s)\n", $1, $3, len[$1] > "/dev/stderr"; bad = 1; exit 1 }
            if (nuse && !($1 in use)) next
            print $1, $2, $3
        }
        END { if (bad) exit 1 }
    ' "$FAI" "$CALLABLE_BED" > "$REGIONS.unmerged"
    REGION_SOURCE="${REGION_SOURCE:+$REGION_SOURCE, }callable BED $(readlink -f "$CALLABLE_BED")"
else
    awk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { len[$1] = $2; next } { print $1, 0, len[$1] }' "$FAI" "$ALLOWED_CONTIGS" > "$REGIONS.unmerged"
fi
# Sort and merge overlapping/adjacent regions.
sort -k1,1 -k2,2n "$REGIONS.unmerged" | awk -F'\t' -v OFS='\t' '
    $1 != c || $2 > e { if (c != "") print c, s, e; c = $1; s = $2; e = $3; next }
    $3 > e { e = $3 }
    END { if (c != "") print c, s, e }
' > "$REGIONS"
rm -f "$REGIONS.unmerged"
[[ -s "$REGIONS" ]] || { echo "ERROR: no callable regions left after applying -c / -b" >&2; exit 1; }
cut -f1 "$REGIONS" | awk '!seen[$0]++' > "$ALLOWED_CONTIGS"
echo "Consensus restricted to: $REGION_SOURCE"
echo "  $(wc -l < "$ALLOWED_CONTIGS") contig(s), $(wc -l < "$REGIONS") region(s), $(awk -F'\t' '{ n += $3 - $2 } END { printf "%.0f", n }' "$REGIONS") bp"
: > "$INT_DIR/outside_regions.tsv"

######################
## Helper functions ##
######################

# Left-align indels and split multiallelic sites; plain-text output.
normalize() {
    bcftools norm -f "$REF_FASTA" -m -any -Ov -o "$2" "$1"
}

# Split a VCF into SNV and indel files. Multiallelic records go to $4.
# Equal-length multi-base substitutions (MNVs) are copied unchanged to $5 and
# decomposed into single-base SNVs (positions where REF == ALT are skipped).
# A derived SNV keeps the caller's FILTER, but its QUAL and INFO are replaced by
# MNV_DERIVED;MNV_EVENT=<POS>:<REF>><ALT>: the caller computed those scores for
# the whole MNV, not for the single base, so they must not be reused per SNV.
split_variants() {
    local in=$1 snv_out=$2 indel_out=$3 multi_out=$4 mnv_out=$5
    : > "$snv_out"; : > "$indel_out"; : > "$multi_out"; : > "$mnv_out"
    awk -F'\t' -v OFS='\t' -v snv="$snv_out" -v indel="$indel_out" -v multi="$multi_out" -v mnv="$mnv_out" '
        /^#CHROM/ {
            print "##INFO=<ID=MNV_DERIVED,Number=0,Type=Flag,Description=\"SNV decomposed from an MNV call; the caller'"'"'s event-level QUAL and INFO were removed\">" > snv
            print "##INFO=<ID=MNV_EVENT,Number=1,Type=String,Description=\"Original MNV as POS:REF>ALT\">" > snv
        }
        /^#/                      { print > snv; print > indel; print > mnv; next }
        $5 ~ /,/                  { print > multi; next }
        length($4) != length($5)  { print > indel; next }
        length($4) == 1           { print > snv; next }
        {
            print > mnv
            ref = $4; alt = $5; pos = $2
            event = pos ":" ref ">" alt
            for (i = 1; i <= length(ref); i++) {
                r = substr(ref, i, 1); a = substr(alt, i, 1)
                if (r == a) continue
                $2 = pos + i - 1; $4 = r; $5 = a; $6 = "."
                $8 = "MNV_DERIVED;MNV_EVENT=" event
                print > snv
            }
        }
    ' "$in"
}

# Print one caller's records in long format: CHROM POS REF ALT TAG FILTER INFO QUAL.
# Only records whose POS lies inside the callable regions are printed; the others
# are counted per contig in $INT_DIR/outside_regions.tsv.
# A FILTER of "." (filters not applied) is recorded as "MISSING" so it is never
# mistaken for "caller absent" or PASS.
emit_long() {
    local tag=$1 file=$2 vtype=$3
    awk -F'\t' -v OFS='\t' -v tag="$tag" -v vtype="$vtype" -v outside="$INT_DIR/outside_regions.tsv" '
        # Regions are sorted and merged; binary search for the last start < POS.
        function inside(c, pos,    lo, hi, mid) {
            if (!(c in NR_C)) return 0
            lo = 1; hi = NR_C[c]
            while (lo < hi) { mid = int((lo + hi + 1) / 2); if (RS_[c, mid] < pos) lo = mid; else hi = mid - 1 }
            return (RS_[c, lo] < pos && pos <= RE_[c, lo])
        }
        FILENAME == ARGV[1] { n = ++NR_C[$1]; RS_[$1, n] = $2 + 0; RE_[$1, n] = $3 + 0; next }
        /^#/ || $5 == "." { next }
        !inside($1, $2 + 0) { out[$1]++; next }
        { print $1, $2, $4, $5, tag, ($7 == "." ? "MISSING" : $7), $8, $6 }
        END { for (c in out) print tag, vtype, c, out[c] >> outside }
    ' "$REGIONS" "$file"
}

# Build the per-variant table from long-format records on stdin.
# Columns: 1 CHROM 2 POS 3 ID 4 REF 5 ALT 6 DECISION
#          7-10  FILTER for mt sa st ms ("." = caller did not report it)
#          11-14 INFO   for mt sa st ms
#          15    CALLER_RECORDS: callers with more than one record, e.g. mt:2 ("." if none)
# One caller counts once per variant. When it has several records, one is kept by
# a fixed order that does not depend on input order (see better() below).
aggregate() {
    awk -F'\t' -v OFS='\t' '
        BEGIN { ntag = split("mt sa st ms", T, " "); NOSCORE = -1e300 }
        # Caller confidence: Mutect2 TLOD, Strelka2 SomaticEVS, otherwise QUAL.
        function score(tag, info, qual,    key, n, i, kv, p) {
            key = (tag == "mt" ? "TLOD" : (tag == "st" ? "SomaticEVS" : ""))
            if (key != "") {
                n = split(info, kv, ";")
                for (i = 1; i <= n; i++) { split(kv[i], p, "="); if (p[1] == key && p[2] != "") return p[2] + 0 }
                return NOSCORE
            }
            return (qual ~ /^[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ ? qual + 0 : NOSCORE)
        }
        function derived(info) { return (info ~ /(^|;)MNV_DERIVED(;|$)/) }
        # 1 if record A should be kept over record B:
        # PASS first, then directly called over MNV-derived, then higher score,
        # then the lexically smaller FILTER/INFO/QUAL text (a total order).
        function better(fa, ia, qa, sa, fb, ib, qb, sb) {
            if ((fa == "PASS") != (fb == "PASS")) return (fa == "PASS")
            if (derived(ia) != derived(ib)) return !derived(ia)
            if (sa != sb) return (sa > sb)
            return ((fa "\t" ia "\t" qa) < (fb "\t" ib "\t" qb))
        }
        {
            k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            if (!(k in seen)) { seen[k] = 1; keys[++nk] = k }
            c = k SUBSEP $5
            NREC[c]++
            sc = score($5, $7, $8)
            if (!(c in F) || better($6, $7, $8, sc, F[c], I[c], Q[c], SC[c])) { F[c] = $6; I[c] = $7; Q[c] = $8; SC[c] = sc }
        }
        END {
            for (n = 1; n <= nk; n++) {
                split(keys[n], f, SUBSEP)
                line = f[1] OFS f[2] OFS "." OFS f[3] OFS f[4] OFS "."
                dup = ""
                for (j = 1; j <= ntag; j++) { c = keys[n] SUBSEP T[j]; line = line OFS ((c in F) ? F[c] : ".") }
                for (j = 1; j <= ntag; j++) { c = keys[n] SUBSEP T[j]; line = line OFS ((c in I) ? I[c] : ".")
                                              if (NREC[c] > 1) dup = dup (dup == "" ? "" : ",") T[j] ":" NREC[c] }
                print line, (dup == "" ? "." : dup)
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
# MNV-derived SNV records carry no TLOD/SomaticEVS (those were computed for the
# whole MNV), so a caller supporting a variant only through a decomposed MNV is
# neither checked nor demoted here.
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
        printf '#CHROM\tPOS\tID\tREF\tALT\tSCENARIO%s_DECISION\tMUTECT2_FILTER\tSAGE_FILTER\tSTRELKA_FILTER\tMUSE_FILTER\tMUTECT2_INFO\tSAGE_INFO\tSTRELKA_INFO\tMUSE_INFO\tCALLER_RECORDS\n' "$sc"
        cat "$table"
    } > "$tsv"

    awk -F'\t' -v OFS='\t' -v sc="$sc" -v pon="$PON_VETO" -v ref="$(readlink -f "$REF_FASTA")" -v fai="$FAI" -v allowed="$ALLOWED_CONTIGS" -v regions="$REGION_SOURCE" '
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
            cl = ""; pc = ""; cf = ""; extra = ""; md = ""; n = 0; p = 0
            for (j = 1; j <= 4; j++) {
                f = $(6 + j)
                if (f == ".") continue
                n++; cl = cl (cl == "" ? "" : ",") T[j]
                g = f; gsub(";", "|", g); cf = cf (cf == "" ? "" : ",") T[j] ":" g
                if (f == "PASS") { p++; pc = pc (pc == "" ? "" : ",") T[j] }
                info = $(10 + j)
                if (info ~ /(^|;)MNV_DERIVED(;|$)/) md = md (md == "" ? "" : ",") T[j]
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
            out = out ";caller_filters=" cf
            if (md != "") out = out ";MNV_DERIVED;mnv_derived_callers=" md
            if ($15 != "" && $15 != ".") out = out ";caller_records=" $15
            out = out extra
            filt = (p >= 2 ? "PASS" : (p == 1 ? "OnePassCaller" : "NoPassCaller"))
            rec[++nr] = $1 OFS $2 OFS "." OFS $4 OFS $5 OFS "." OFS filt OFS out
        }
        END {
            print "##fileformat=VCFv4.2"
            print "##source=consensus_scenarios.sh scenario " sc
            print "##reference=file://" ref
            print "##consensus_regions=" regions
            # ##contig lines for the allowed contigs only, with lengths from the reference index.
            while ((getline line < allowed) > 0) use[line] = 1
            while ((getline line < fai) > 0) { split(line, fl, "\t"); if (fl[1] in use) print "##contig=<ID=" fl[1] ",length=" fl[2] ">" }
            print "##FILTER=<ID=PASS,Description=\"Called PASS by >= 2 callers" (sc == 6 ? "; also passed scenario 6 bias filter and TLOD/SomaticEVS checks" (pon == 1 ? " and the Mutect2 panel-of-normals veto" : "") : "") "\">"
            print "##FILTER=<ID=OnePassCaller,Description=\"Called PASS by exactly 1 caller (see caller_filters)\">"
            print "##FILTER=<ID=NoPassCaller,Description=\"Not called PASS by any caller (see caller_filters)\">"
            print "##INFO=<ID=callers,Number=.,Type=String,Description=\"Callers reporting the variant, any FILTER (mt=Mutect2, sa=SAGE, st=Strelka2, ms=MuSE)\">"
            print "##INFO=<ID=ncallers,Number=1,Type=Integer,Description=\"Number of distinct callers reporting the variant\">"
            print "##INFO=<ID=pass_callers,Number=.,Type=String,Description=\"Callers with FILTER=PASS (scenario 6: after threshold demotion)\">"
            print "##INFO=<ID=npass,Number=1,Type=Integer,Description=\"Number of distinct PASS callers\">"
            print "##INFO=<ID=caller_filters,Number=.,Type=String,Description=\"Per-caller FILTER as caller:filter; multiple filters joined with |\">"
            print "##INFO=<ID=MNV_DERIVED,Number=0,Type=Flag,Description=\"At least one caller supports this SNV only through a decomposed MNV call\">"
            print "##INFO=<ID=mnv_derived_callers,Number=.,Type=String,Description=\"Callers whose record is an SNV decomposed from their MNV call (event-level scores removed)\">"
            print "##INFO=<ID=caller_records,Number=.,Type=String,Description=\"Callers with more than one record for this variant, as caller:count; the kept record follows PASS > direct call > higher score > lexical order\">"
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
    "$INT_DIR/${S}_mutect_snv.vcf" "$INT_DIR/${S}_mutect_indel.vcf" "$SPECIAL_DIR/${S}_mutect.multiMutation.txt" \
    "$SPECIAL_DIR/${S}_mutect.mnv.vcf"

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
    "$INT_DIR/${S}_sage_snv.vcf" "$INT_DIR/${S}_sage_indel.vcf" "$SPECIAL_DIR/${S}_sage.multiMutation.txt" \
    "$SPECIAL_DIR/${S}_sage.mnv.vcf"

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
    emit_long mt "$INT_DIR/${S}_mutect_snv.vcf" snv
    emit_long sa "$INT_DIR/${S}_sage_snv.vcf" snv
    emit_long st "$INT_DIR/${S}_strelka_snv.vcf" snv
    emit_long ms "$INT_DIR/${S}_muse_snv.vcf" snv
} | aggregate > "$INT_DIR/${S}_snv_table.tsv"

{
    emit_long mt "$INT_DIR/${S}_mutect_indel.vcf" indel
    emit_long sa "$INT_DIR/${S}_sage_indel.vcf" indel
    emit_long st "$INT_DIR/${S}_strelka_indel.vcf" indel
} | aggregate > "$INT_DIR/${S}_indel_table.tsv"

echo "  Distinct SNVs: $(wc -l < "$INT_DIR/${S}_snv_table.tsv"), distinct indels: $(wc -l < "$INT_DIR/${S}_indel_table.tsv")"

# Records outside the callable regions never enter the consensus; list them per caller and contig.
{
    printf 'caller\tvariant_type\tcontig\trecords\n'
    sort -t$'\t' -k3,3V -k1,1 -k2,2 "$INT_DIR/outside_regions.tsv"
} > "$SPECIAL_DIR/${S}_outside_callable_regions.tsv"
n_out=$(awk -F'\t' '{ n += $4 } END { print n + 0 }' "$INT_DIR/outside_regions.tsv")
echo "  Records outside the callable regions (excluded): $n_out, on $(cut -f3 "$INT_DIR/outside_regions.tsv" | sort -u | wc -l) contig(s); see special_mutations/${S}_outside_callable_regions.tsv"

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
