#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_bam> <assembly> <min_cov> <min_len> <output_folder>"
    exit 1
fi

# Assign input arguments to variables
input_bam="$1"
assembly="$2"
min_cov="$3"
min_len="$4"
output_folder="$5"

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Extract the file name without the extension
filename=$(basename "${input_bam%.bam}")
filename=$(basename "${filename%.sorted}")

# Run samtools depth
samtools depth -a "$input_bam" > "$output_folder/${filename}.bedgraph"

# Filter for low coverage
awk '$3 < '$min_cov "$output_folder/${filename}.bedgraph" > "$output_folder/${filename}-low_coverage.bedgraph"

# Create low coverage TSV
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $2}' "$output_folder/${filename}-low_coverage.bedgraph" > "$output_folder/${filename}-low_coverage.bed"

# Merge low coverage intervals
bedtools merge -i "$output_folder/${filename}-low_coverage.bed" -d 100 > "$output_folder/${filename}-low_coverage_merged.tmp.bed"
# Remove one base from the start and end of each interval
awk 'BEGIN{OFS="\t"}{$2=$2+1; $3=$3-1}1' "$output_folder/${filename}-low_coverage_merged.tmp.bed" > "$output_folder/${filename}-low_coverage_merged.bed"
# Remove the temporary file
rm "$output_folder/${filename}-low_coverage_merged.tmp.bed"

# Remove small sequences
awk '($3 - $2) >= '$min_len "$output_folder/${filename}-low_coverage_merged.bed" > "$output_folder/${filename}-GD.bed"

# Extract FASTA sequences
bedtools getfasta -fi "$assembly" -bed "$output_folder/${filename}-GD.bed" -fo "$output_folder/${filename}-GD.fasta"

rm "$output_folder/${filename}.bedgraph"
rm "$output_folder/${filename}-low_coverage.bedgraph"
rm "$output_folder/${filename}-low_coverage.bed"
rm "$output_folder/${filename}-low_coverage_merged.bed"
