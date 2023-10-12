#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 4 ]; then
  echo "Usage: $0 <fastq_file> <assembly> <mapped_folder> <threads>"
  exit 1
fi

# Assign the input arguments to variables
fastq=$1
assembly=$2
mapped_folder=$3
thr=$4

# Check if the folder exists
if [ ! -f "$fastq" ]; then
  echo "Fastq file does not exist: $fastq_file"
  exit 1
fi

# Check if the fasta file exists
if [ ! -f "$assembly" ]; then
  echo "Assembly file does not exist: $assembly"
  exit 1
fi

# Extract the file name without the extension
filename=$(basename "$fastq" .fastq.gz)

# Run BWA-MEM to map fastq to fasta
bwa mem -t "${thr}" "${assembly}" "${fastq}" > "${mapped_folder}/${filename}.sam"
samtools view -bS -o "${mapped_folder}/${filename}.bam" "${mapped_folder}/${filename}.sam" > "/dev/null"
samtools sort "${mapped_folder}/${filename}.bam" -o "${mapped_folder}/${filename}.sorted.bam"
samtools index "${mapped_folder}/${filename}.sorted.bam"
rm "${mapped_folder}/${filename}.sam"
rm "${mapped_folder}/${filename}.bam"
done

echo "${filename} mapped successfully to ${assembly}"
echo "Extracting low coverage sequences from ${filename}"
bash "scripts/bam2fasta.sh" "${mapped_folder}/${filename}.sorted.bam" "${assembly}" "${mapped_folder}"

blastn -query "${mapped_folder}/unmapped.fasta" -subject "${mapped_folder}/unmapped.fasta" -out "${mapped_folder}/unmapped.blast"  -outfmt 6
awk '($12) >= 1000' "${mapped_folder}/unmapped.blast" > "${mapped_folder}/unmapped-filtered.blast"
mkdir "${mapped_folder}/clusters/"
echo "Finding repetitive sequences..."
python "scripts/blast2clusters.py" "${mapped_folder}/unmapped-filtered.blast" "${mapped_folder}/unmapped.fasta" "${mapped_folder}/clusters/"

echo "Extracting consensus sequences of the invaders..."
for fasta in "${mapped_folder}/clusters/"*.fasta; do
    # Define the output file name based on the input file name
    output_MSA="${fasta%.fasta}.MSA"
    output_consensus="${fasta%.fasta}.consensus"
    output_standard="${fasta%.fasta}"
    
    # Refine clusters
    #python "scripts/refine-clusters.py" "$fasta" "$output_standard" "$assembly"

    # Run MUSCLE with input and output files
    muscle -in "${output_standard}.fasta" -out "$output_MSA"
    python "scripts/MSA2consensus.py" "$output_MSA" "$output_consensus"
done

cat "${mapped_folder}/clusters/"*consensus > "${mapped_folder}/candidates.fasta"
echo "Extracting consensus sequences of the invaders..."