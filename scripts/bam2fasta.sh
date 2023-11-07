#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_bam> <assembly> <min_cov> <min_len> <max_d> <output_folder>"
    exit 1
fi

# Assign input arguments to variables
input_bam="$1"
assembly="$2"
min_cov="$3"
min_len="$4"
max_d="$5"
output_folder="$6"

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
bedtools merge -i "$output_folder/${filename}-low_coverage.bed" -d "$max_d" > "$output_folder/${filename}-low_coverage_merged.tmp.bed"
# Remove one base from the start and end of each interval
awk 'BEGIN{OFS="\t"}{$2=$2+1; $3=$3-1}1' "$output_folder/${filename}-low_coverage_merged.tmp.bed" > "$output_folder/${filename}-low_coverage_merged.bed"
# Remove the temporary file
rm "$output_folder/${filename}-low_coverage_merged.tmp.bed"

# Remove small sequences
awk '($3 - $2) >= '$min_len "$output_folder/${filename}-low_coverage_merged.bed" > "$output_folder/${filename}-GD.bed"

# Assign credibility to each sequence based on genomic context coverage
current_dir=$(dirname "$(readlink -f "$0")")
bash "$current_dir/credibility.sh" "$output_folder/${filename}-GD.bed" "$output_folder/${filename}.bedgraph" "$input_bam" "$assembly" "$output_folder"
python "$current_dir/credibility.py" "$output_folder/${filename}-GD.bed" "$output_folder/${filename}-GD-flanking.credibility" "$output_folder/${filename}-GD.fai"

# Use a loop to process each line in the BED file and get the corresponding cred_value
bedtools getfasta -fi "$assembly" -bed "$output_folder/${filename}-GD-credibility.bed" -fo "$output_folder/${filename}-GD.fasta" -name
sed -i'' -e 's/::.*$//' "$output_folder/${filename}-GD.fasta"  

# Remove the temporary files
rm "$output_folder/${filename}.bedgraph"
rm "$output_folder/${filename}-low_coverage.bedgraph"
rm "$output_folder/${filename}-low_coverage.bed"
rm "$output_folder/${filename}-low_coverage_merged.bed"