#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <bed> <bedgraph> <bam> <fasta_genome> <output_folder>"
    exit 1
fi

# Assign input arguments to variables
bed="$1"
bedgraph="$2"
bam="$3"
fasta_genome="$4"
output_folder="$5"

# Extract the file name without the extension
filename=$(basename "${bed%.bed}")

# Function to generate genome file from FASTA assembly
samtools faidx "$fasta_genome" && mv "$fasta_genome.fai" "$output_folder/${filename}.fai"

# Slop 10,000 nucleotides before the regions
bedtools flank -i "$bed" -g "$output_folder/${filename}.fai" -b 10000 |
awk -F'\t' -v OFS='\t' '{$2 = ($2 == 0 ? 1 : $2); print}' |
bedtools sort -i - > "$output_folder/${filename}-flanking.bed"

bedtools coverage -sorted -g "$output_folder/${filename}.fai" -a "$output_folder/${filename}-flanking.bed" -b "$bam" -mean > "$output_folder/${filename}-flanking.bedgraph"

# Calculate mean coverage of the whole genome
mean_coverage_whole_genome=$(awk '{ total += $3; count++ } END { if(count > 0) print total/count }' "$bedgraph")

awk -v mean="$mean_coverage_whole_genome" '{ $5 = $4 / mean; print }' "$output_folder/${filename}-flanking.bedgraph" > "$output_folder/${filename}-flanking.credibility"