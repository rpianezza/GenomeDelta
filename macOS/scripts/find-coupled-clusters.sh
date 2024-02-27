#!/bin/bash

# Assign input arguments to variables
fasta_folder="$1"
output_folder="$2"
max_d="$3"
assembly="$4"
input_bam="$5"
bedgraph="$6"
prefix="$7"

# Get the directory of the currently running shell script
script_dir="$(dirname "$0")"

for fasta_file in "$fasta_folder"/*.fasta; do
    # Extract the filename without extension
    filename="$prefix"

    # Remove the last part after the last underscore
    bed_file_name="${filename}"

    # Set the bed file path
    bed_file="$output_folder/$bed_file_name.bed"

    # Run the Python script with the arguments
    python "$script_dir/find-coupled-clusters.py" "$fasta_file" "$max_d" "$bed_file"
done

# Create an associative array to store file hashes
declare -afFirtx file_hashes

# Function to calculate MD5 hash of a file
calculate_hash() {
    md5 "$1" | awk '{print $1}'
}

# Iterate over all files in the folder
for file_path in "$output_folder"/*txt; do
    # Check if the file is a regular file
    if [ -f "$file_path" ]; then
        # Calculate the hash of the file content
        file_hash=$(calculate_hash "$file_path")

        # Check if the hash is already in the array
        if [ "${file_hashes[$file_hash]}" ]; then
            # If a file with the same hash exists, remove the redundant copy
            echo "Removing redundant copy: $file_path"
            rm "$file_path"
        else
            # Otherwise, store the hash in the array
            file_hashes["$file_hash"]=$file_path
        fi
    fi
done

for txt_file in "$output_folder"/*.txt; do
    # Extract the filename without extension
    filename=$(basename -- "$txt_file")
    filename_no_extension="${filename%.*}"

    # Set the bed file path
    bed_file="$output_folder/$filename_no_extension.bed"

    # Convert the txt file to a BED file
    awk -F '[:-]' '{print substr($1, 2) "\t" $2 "\t" $3}' "$txt_file" > "$bed_file"

    # Assign credibility to each sequence based on genomic context coverage
    bash "$script_dir/credibility.sh" "$bed_file" "$bedgraph" "$input_bam" "$assembly" "$output_folder"
    python "$script_dir/credibility.py" "$bed_file" "$output_folder/${filename_no_extension}-flanking.credibility" "$output_folder/${filename_no_extension}.fai"
    
    # Run bedtools getfasta for each BED file
    fasta_output="$output_folder/$filename_no_extension.fasta"
    bedtools getfasta -fi "$assembly" -bed "$output_folder/${filename_no_extension}-credibility.bed" -fo "$fasta_output" -name
    sed -i'' -e 's/::.*$//' "$fasta_output"

    muscle -in "$fasta_output" -out "$output_folder/$filename_no_extension.MSA"
    python "$script_dir/MSA2consensus.py" "$output_folder/$filename_no_extension.MSA" "$output_folder/$filename_no_extension.consensus"
    rm "$output_folder/$filename_no_extension.MSA"
    rm "$output_folder/$filename_no_extension.fasta-e"
    rm "$output_folder/$filename_no_extension-flanking.bed"
    rm "$output_folder/$filename_no_extension-credibility.bed"
    rm "$output_folder/$filename_no_extension-flanking.credibility"
    rm "$output_folder/$filename_no_extension-flanking.bedgraph"
    rm "$output_folder/$filename_no_extension.bed"
    rm "$output_folder/$filename_no_extension.txt"
    rm "$output_folder/$filename_no_extension.fai"
done