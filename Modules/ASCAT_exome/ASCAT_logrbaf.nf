nextflow.enable.dsl=2

process ASCAT_logrbaf {
    conda "${params.ascat_env}"
    label 'process_low'
    publishDir("${params.ascat_dir}", mode: 'copy', overwrite: true, failOnError: false)

    errorStrategy 'retry'
    maxRetries 3

    input:
    // Accept separate collections for normal and tumor files
    tuple val(patient), val(sample), val(gender), path(normal_files), path(tumor_files)

    output:
    tuple val(patient),
          val(sample),
          val(gender),
          path("${patient}_${sample}_tumor_tumourLogR.txt"),
          path("${patient}_${sample}_tumor_tumourBAF.txt"),
          path("${patient}_${sample}_tumor_normalLogR.txt"),
          path("${patient}_${sample}_tumor_normalBAF.txt"),
          emit: ascat_input

    script:
    def min_depth = params.min_depth ?: 20  // Default to 20 if not defined in params

    """
    #!/usr/bin/env python3
    import os
    import glob
    import math

    # Output files
    tumor_logr_file = "${patient}_${sample}_tumor_tumourLogR.txt"
    tumor_baf_file = "${patient}_${sample}_tumor_tumourBAF.txt"
    normal_logr_file = "${patient}_${sample}_tumor_normalLogR.txt"
    normal_baf_file = "${patient}_${sample}_tumor_normalBAF.txt"

    # Initialize files with headers
    with open(tumor_logr_file, 'w') as f1:
        f1.write("ID\\tChromosome\\tPosition\\t${patient}_${sample}_tumor\\n")
    with open(tumor_baf_file, 'w') as f2:
        f2.write("ID\\tChromosome\\tPosition\\t${patient}_${sample}_tumor\\n")
    with open(normal_logr_file, 'w') as f3:
        f3.write("ID\\tChromosome\\tPosition\\t${patient}_${sample}_normal\\n")
    with open(normal_baf_file, 'w') as f4:
        f4.write("ID\\tChromosome\\tPosition\\t${patient}_${sample}_normal\\n")

    # Get lists of normal and tumor files
    normal_files = sorted(glob.glob("*normal*"))
    tumor_files = sorted(glob.glob("*tumor*"))

    # Organize files by chromosome
    chr_files = {}
    for nf in normal_files:
        chr_num = os.path.basename(nf).split('_')[-1].replace('.txt', '')
        if chr_num not in chr_files:
            chr_files[chr_num] = {'normal': None, 'tumor': None}
        chr_files[chr_num]['normal'] = nf

    for tf in tumor_files:
        chr_num = os.path.basename(tf).split('_')[-1].replace('.txt', '')
        if chr_num not in chr_files:
            chr_files[chr_num] = {'normal': None, 'tumor': None}
        chr_files[chr_num]['tumor'] = tf

    # Batch processing of positions to improve writing efficiency
    batch_size = 1000
    total_positions = 0

    # Define a function to write a batch of results to files
    def write_batch_to_files(tumor_logr_batch, tumor_baf_batch, normal_logr_batch, normal_baf_batch):
        with open(tumor_logr_file, 'a') as f1:
            f1.write(''.join(tumor_logr_batch))
        with open(tumor_baf_file, 'a') as f2:
            f2.write(''.join(tumor_baf_batch))
        with open(normal_logr_file, 'a') as f3:
            f3.write(''.join(normal_logr_batch))
        with open(normal_baf_file, 'a') as f4:
            f4.write(''.join(normal_baf_batch))

    # Process each chromosome
    for chr_num, files in chr_files.items():
        normal_file = files.get('normal')
        tumor_file = files.get('tumor')

        if not normal_file or not tumor_file:
            continue

        # Read normal file into memory
        normal_data = {}
        with open(normal_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('CHR') or line.startswith('chrom'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) >= 7:  # Make sure we have enough columns
                    chrom, pos, a, c, g, t, depth = parts[0:7]

                    # Convert to integers
                    try:
                        a, c, g, t, depth = int(a), int(c), int(g), int(t), int(depth)
                    except ValueError:
                        continue

                    # Filter by minimum depth
                    if depth >= ${min_depth}:
                        normal_data[f'{chrom}_{pos}'] = (a, c, g, t, depth)

        # Process tumor file and generate output
        tumor_logr_batch = []
        tumor_baf_batch = []
        normal_logr_batch = []
        normal_baf_batch = []

        with open(tumor_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('CHR') or line.startswith('chrom'):
                    continue
                parts = line.strip().split('\\t')
                if len(parts) >= 7:  # Make sure we have enough columns
                    chrom, pos, a, c, g, t, depth = parts[0:7]

                    # Convert to integers
                    try:
                        a, c, g, t, depth = int(a), int(c), int(g), int(t), int(depth)
                    except ValueError:
                        continue

                    # Filter by minimum depth
                    if depth < ${min_depth}:
                        continue

                    key = f'{chrom}_{pos}'

                    if key in normal_data:
                        # Get counts from normal data
                        n_a, n_c, n_g, n_t, n_depth = normal_data[key]

                        # Calculate totals
                        tumor_total = a + c + g + t
                        normal_total = n_a + n_c + n_g + n_t

                        # Skip positions with no coverage
                        if tumor_total == 0 or normal_total == 0:
                            continue

                        # Determine reference and alternate alleles
                        # Find the most frequent base in normal sample
                        normal_bases = [n_a, n_c, n_g, n_t]
                        normal_max_idx = normal_bases.index(max(normal_bases))
                        bases = ['a', 'c', 'g', 't']
                        ref_base = bases[normal_max_idx]

                        # Get reference counts based on ref_base
                        if ref_base == 'a':
                            tumor_ref_count = a
                            normal_ref_count = n_a
                        elif ref_base == 'c':
                            tumor_ref_count = c
                            normal_ref_count = n_c
                        elif ref_base == 'g':
                            tumor_ref_count = g
                            normal_ref_count = n_g
                        elif ref_base == 't':
                            tumor_ref_count = t
                            normal_ref_count = n_t

                        # Calculate alternate counts
                        tumor_alt_count = tumor_total - tumor_ref_count
                        normal_alt_count = normal_total - normal_ref_count

                        # Clean chromosome name (remove 'chr' prefix if present)
                        chr_clean = chrom.replace('chr', '')

                        # Create unique ID for this position
                        pos_id = f'{chr_clean}_{pos}'

                        # Calculate BAF and LogR values
                        tumor_baf_value = tumor_alt_count / tumor_total if tumor_total > 0 else 0
                        normal_baf_value = normal_alt_count / normal_total if normal_total > 0 else 0
                        logr_value = math.log2(tumor_total / normal_total) if normal_total > 0 else 0

                        # Add to batches with 4-column format (ID, Chromosome, Position, Value)
                        tumor_logr_batch.append(f'{pos_id}\\t{chr_clean}\\t{pos}\\t{logr_value}\\n')
                        tumor_baf_batch.append(f'{pos_id}\\t{chr_clean}\\t{pos}\\t{tumor_baf_value}\\n')
                        normal_logr_batch.append(f'{pos_id}\\t{chr_clean}\\t{pos}\\t0\\n')
                        normal_baf_batch.append(f'{pos_id}\\t{chr_clean}\\t{pos}\\t{normal_baf_value}\\n')

                        total_positions += 1

                        # Write batch when it reaches the specified size
                        if len(tumor_logr_batch) >= batch_size:
                            write_batch_to_files(tumor_logr_batch, tumor_baf_batch, normal_logr_batch, normal_baf_batch)
                            tumor_logr_batch = []
                            tumor_baf_batch = []
                            normal_logr_batch = []
                            normal_baf_batch = []

        # Write any remaining data in batches
        if tumor_logr_batch:
            write_batch_to_files(tumor_logr_batch, tumor_baf_batch, normal_logr_batch, normal_baf_batch)
        
        # Clear normal_data to free memory
        normal_data.clear()

    print(f"Completed processing for sample ${patient}_${sample}")
    """
}
