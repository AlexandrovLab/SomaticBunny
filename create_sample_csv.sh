#!/bin/bash

# Define the path to the CSV output
output_file="sample.csv"

# Write the header to the output CSV
echo "patient,status,fastq_1,fastq_2,type,file" > "$output_file"

# Define the base path for the fastq files
base_path="/your/fastq/path/"

# Loop over each name in list.txt
for patient in $(cat /tscc/lustre/restricted/alexandrov-ddn/users/tiy002/dbGaP/NB_TARGET/ecDNA_analysis/BWA/list.txt); do
    # Add normal and tumor samples for each patient
    echo "$patient,normal,$base_path/"$patient"_normal_1.fastq.gz,$base_path/"$patient"_normal_2.fastq.gz,genome,fastq" >> "$output_file"
    echo "$patient,tumor,$base_path/"$patient"_tumor_1.fastq.gz,$base_path/"$patient"_tumor_2.fastq.gz,genome,fastq" >> "$output_file"
done

echo "sample.csv created successfully."
