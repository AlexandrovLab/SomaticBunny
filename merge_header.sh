set -eu

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 first.vcf [other.vcf ...] > merged_header.txt" >&2
    exit 1
fi

for f in "$@"; do
    [[ -r "$f" ]] || { echo "ERROR: cannot read $f" >&2; exit 1; }
done

echo "##fileformat=VCFv4.2"

# Diagnostic: how many INFO lines each input file contributes (stderr)
for f in "$@"; do
    n=$(zcat -f "$f" | awk '/^#CHROM/ {exit} !/^#/ {exit} /^##INFO=/ {c++} END {print c+0}')
    echo "  $f: $n INFO lines" >&2
done

# INFO: merged across all files
for f in "$@"; do
    zcat -f "$f" | awk -v src="$f" '
        /^#CHROM/ { exit }
        !/^#/     { exit }
        /^##INFO=/ { print src "\t" $0 }
    '
done | awk -F'\t' '
{
    src = $1; line = $2
    if (match(line, /ID=[^,>]+/)) id = substr(line, RSTART + 3, RLENGTH - 3)
    else { print "WARNING: no ID in line from " src ": " line > "/dev/stderr"; next }

    if (!(id in seen)) {
        seen[id] = line; from[id] = src; order[++n] = id
    } else if (seen[id] != line) {
        printf "WARNING: INFO/%s differs in %s (kept definition from %s)\n", id, src, from[id] > "/dev/stderr"
    }
}
END {
    for (i = 1; i <= n; i++) print seen[order[i]]
    printf "Merged %d unique INFO entries\n", n > "/dev/stderr"
}'

echo "##INFO=<ID=callers,Number=.,Type=String,Description=\"Variant callers that called the variant\">"
echo "##INFO=<ID=caller_filtering_status,Number=.,Type=String,Description=\"Filtering status reported by each caller\">"


# contig: from the first file only
zcat -f "$1" | awk '
    /^#CHROM/ { exit }
    !/^#/     { exit }
    /^##contig=/ { print }
'
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
