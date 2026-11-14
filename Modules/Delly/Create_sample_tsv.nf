process Create_sample_tsv {
    publishDir("${params.Delly_dir}", mode: 'copy')
    
    input:
    path(sample_csv)
    
    output:
    path("samples.tsv"), emit: samples_tsv
    
    script:
    """
    #!/usr/bin/env python3
    
    import csv
    
    # Read the sample sheet and create the samples.tsv
    with open("${sample_csv}", 'r') as f, open("samples.tsv", 'w') as out:
        reader = csv.DictReader(f)
        for row in reader:
            patient = row['patient']
            sample = row['sample']
            status = row['status']
            
            # Create a sample ID as it would appear in the VCF/BCF according to Delly
            sample_id = f"{patient}_{sample}_{status}"
            
            # Map status to Delly's terminology
            sample_type = "tumor" if status == "tumor" else "control"
            
            # Write to samples.tsv
            out.write(f"{sample_id}\\t{sample_type}\\n")
    """
}
