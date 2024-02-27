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

# Slop 10,000 nucleotides before the regions
bedtools flank -i "$output_folder/${filename}-GD.bed" -g "$output_folder/${filename}.fai" -b 10000 |
awk -F'\t' -v OFS='\t' '{$2 = ($2 == 0 ? 1 : $2); print}' |
bedtools sort -i - > "$output_folder/${filename}-flanking.bed"
#samtools depth -a -b "$output_folder/${filename}-flanking.bed" "$bam" | awk 'BEGIN{OFS="\t"} {if(prev_contig==$1 && prev_end+1==$2) {sum+=$3; count+=1} else {if(prev_contig!="") print prev_contig, start, prev_end, sum/count; prev_contig=$1; start=$2-1; sum=$3; count=1} prev_end=$2} END{if(prev_contig!="") print prev_contig, start, prev_end, sum/count}' > "$output_folder/${filename}-flanking.bedgraph"

#samtools depth -a -b "$output_folder/${filename}-flanking.bed" "$bam" > "$output_folder/${filename}-flanking.bedgraph.raw"

bedtools coverage -a "$output_folder/${filename}-flanking.bed" -b "$bam" -mean > "$output_folder/${filename}-flanking.bedgraph"

# Calculate mean coverage of the whole genome
mean_coverage_whole_genome=$(awk '{ total += $3; count++ } END { if(count > 0) print total/count }' "$bedgraph")

awk -v mean="$mean_coverage_whole_genome" "{ \$5 = 2 * (\$4 / (\$4 + mean)) - 1; print }" "$output_folder/${filename}-flanking.bedgraph" > "$output_folder/${filename}-GD-flanking.credibility"