#!/usr/bin/env python3
"""
Split a sequence dictionary into N interval lists of WHOLE contigs, for scattering GATK ApplyBQSR.

Every contig of the dictionary (primary chromosomes, alt/decoy/viral/unplaced contigs, chrM) is written to
exactly one list, so a scattered ApplyBQSR run over the N lists emits every mapped and placed-unmapped read
exactly once: no read spans a boundary (boundaries are whole contigs) and no contig is left out. Unplaced
unmapped reads are not covered by any contig and must be requested with "-L unmapped" on one of the chunks.

Groups are balanced by total length (longest-processing-time rule: each contig, largest first, goes to the
currently smallest group), so the largest group is max(longest contig, ~genome/N).

usage: make_contig_groups.py <ref.dict> <N> <out_dir>
writes: <out_dir>/<k>-contigs.interval_list for k = 1..N (Picard interval_list format, header copied from the dict)
"""
import os
import sys


def main(dict_path, n_groups, out_dir):
    header, contigs = [], []
    with open(dict_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.startswith('@'):
                break
            header.append(line)
            if line.startswith('@SQ'):
                fields = dict(tag.split(':', 1) for tag in line.split('\t')[1:] if ':' in tag)
                contigs.append((fields['SN'], int(fields['LN'])))
    if not contigs:
        sys.exit(f"ERROR: no @SQ lines in {dict_path}")
    if n_groups < 1 or n_groups > len(contigs):
        sys.exit(f"ERROR: N must be between 1 and the number of contigs ({len(contigs)}), got {n_groups}")

    order = {name: i for i, (name, _) in enumerate(contigs)}           # dictionary order, kept inside each group
    groups = [[] for _ in range(n_groups)]
    totals = [0] * n_groups
    for name, length in sorted(contigs, key=lambda c: (-c[1], order[c[0]])):
        k = min(range(n_groups), key=lambda i: (totals[i], i))
        groups[k].append(name)
        totals[k] += length

    assigned = [name for g in groups for name in g]
    assert sorted(assigned) == sorted(name for name, _ in contigs), "every contig must be assigned exactly once"

    os.makedirs(out_dir, exist_ok=True)
    lengths = dict(contigs)
    for k, group in enumerate(groups, start=1):
        group.sort(key=lambda name: order[name])
        path = os.path.join(out_dir, f"{k}-contigs.interval_list")
        with open(path, 'w') as out:
            out.write('\n'.join(header) + '\n')
            for name in group:
                out.write(f"{name}\t1\t{lengths[name]}\t+\t.\n")
        print(f"{os.path.basename(path)}\t{len(group)} contig(s)\t{totals[k-1]:,} bp\t{group[0]}{' ...' if len(group) > 1 else ''}")
    print(f"total\t{len(contigs)} contigs\t{sum(lengths.values()):,} bp\tlargest group {max(totals):,} bp, smallest {min(totals):,} bp")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
