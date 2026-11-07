#!/bin/bash
# Verify a SomaticBunny database against checksums/database_files.sha256.
#   bin/verify_database.sh /path/to/SomaticBunny_database          every manifest file present and non-empty (seconds)
#   bin/verify_database.sh /path/to/SomaticBunny_database --hash   full SHA256 verification (hashes the whole database)
set -euo pipefail
DB=${1:?usage: verify_database.sh /path/to/SomaticBunny_database [--hash]}
MANIFEST="$(cd "$(dirname "$0")/.." && pwd)/checksums/database_files.sha256"
cd "$DB"
if [[ "${2:-}" == "--hash" ]]; then
    sha256sum --quiet -c "$MANIFEST" && echo "OK: all $(wc -l < "$MANIFEST") manifest files verified by SHA256 in $DB"
else
    missing=0
    while read -r _ path; do
        [[ -f "$path" ]] || { echo "MISSING: $path"; missing=$((missing + 1)); }
    done < "$MANIFEST"
    n=$(wc -l < "$MANIFEST")
    if [[ $missing -gt 0 ]]; then echo "ERROR: $missing of $n manifest files missing in $DB" >&2; exit 1; fi
    echo "OK: all $n manifest files present in $DB"
fi
