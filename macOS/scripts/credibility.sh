#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <bed> <bedgraph> <bam> <fasta_genome> <output_folder> <prefix>"
    exit 1
fi

# Assign input arguments to variables
bed="$1"
bedgraph="$2"
bam="$3"
fasta_genome="$4"
output_folder="$5"
prefix="$6"

# Extract the file name without the extension
filename="$prefix"

# Function to generate genome file from FASTA assembly
samtools faidx "$fasta_genome" && mv "$fasta_genome.fai" "$output_folder/${filename}.fai"
sort -k1,1 "$output_folder/${filename}.fai" > "$output_folder/${filename}.sorted.fai"
mv "$output_folder/${filename}.sorted.fai" "$output_folder/${filename}.fai"

# Slop 10,000 nucleotides before the regions
bedtools flank -i "$output_folder/${filename}-GD.bed" -g "$output_folder/${filename}.fai" -b 10000 |
awk -F'\t' -v OFS='\t' '{$2 = ($2 == 0 ? 1 : $2); print}' |
bedtools sort -i - > "$output_folder/${filename}-flanking.bed"

sort -k1,1 -k2,2n "$output_folder/${filename}-flanking.bed" > "$output_folder/${filename}-flanking.sorted.bed"
mv "$output_folder/${filename}-flanking.sorted.bed" "$output_folder/${filename}-flanking.bed"
bedtools coverage -sorted -a "$output_folder/${filename}-flanking.bed" -b "$bam" -mean -g "$output_folder/${filename}.fai" > "$output_folder/${filename}-flanking.bedgraph"


# Calculate mean coverage of the whole genome
mean_coverage_whole_genome=$(awk '{ total += $3; count++ } END { if(count > 0) print total/count }' "$bedgraph")

awk -v mean="$mean_coverage_whole_genome" "{ \$5 = 2 * (\$4 / (\$4 + mean)) - 1; print }" "$output_folder/${filename}-flanking.bedgraph" > "$output_folder/${filename}-GD-flanking.credibility"