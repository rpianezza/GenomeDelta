#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_bam> <assembly> <output_folder>"
    exit 1
fi

# Assign input arguments to variables
input_bam="$1"
assembly="$2"
output_folder="$3"

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Run samtools depth
samtools depth -a "$input_bam" > "$output_folder/coverage.bedgraph"

# Filter for low coverage
awk '$3 < 2' "$output_folder/coverage.bedgraph" > "$output_folder/low_coverage.bedgraph"

# Create low coverage TSV
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $2}' "$output_folder/low_coverage.bedgraph" > "$output_folder/low_coverage.bed"

# Merge low coverage intervals
bedtools merge -i "$output_folder/low_coverage.bed" -d 100 > "$output_folder/low_coverage_merged.bed"

# Remove small sequences
awk '($3 - $2) >= 1000' "$output_folder/low_coverage_merged.bed" > "$output_folder/unmapped.bed"

# Extract FASTA sequences
bedtools getfasta -fi "$assembly" -bed "$output_folder/unmapped.bed" -fo "$output_folder/unmapped.fasta"

rm "$output_folder/coverage.bedgraph"
rm "$output_folder/low_coverage.bedgraph"
rm "$output_folder/low_coverage.bed"
rm "$output_folder/low_coverage_merged.bed"
