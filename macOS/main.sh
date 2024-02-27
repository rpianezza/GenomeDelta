#!/bin/bash

# Set default values for optional arguments
d=100 # maximum distance accepted for a gap between 2 low coverage sequences
min_cov=1 # minimum coverage of a position to be considered NON low coverage (in this case, only cov=0 is included in low coverages)
min_len=1000 # minimum length for a low-coverage sequence to be included in the output
min_bitscore=1000 # minimum bitscore in the BLAST alignment to consider the sequences part of the same cluster
refine_d=2500 # maximum distance to check coupled-clusters (e.g. if the option "refine" is activated, the script will merge clusters if their insertions are closer than refine_d)
prefix="file"

# Initialize variables
fastq_set=0
bam_set=0

# Parse command line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fq) 
            fastq="$2"; 
            fastq_set=1; 
            shift 
            ;;
        --bam)
            bam="$2"; 
            bam_set=1; 
            shift 
            ;;
        --fa) 
            assembly="$2"; 
            shift 
            ;;
        --of) 
            mapped_folder="$2"; 
            shift 
            ;;
        --t) 
            thr="$2"; 
            shift 
            ;;
        --d) 
            d="$2"; 
            shift 
            ;;
        --min_cov) 
            min_cov="$2"; 
            shift 
            ;;
        --min_len) 
            min_len="$2"; 
            shift
            ;;
        --min_bitscore) 
            min_bitscore="$2"; 
            shift
            ;;
        --refine) 
            refine_set=1; 
            shift
            ;;
         --refine_d)
            refine_d="$2"; 
            shift
            ;;
         --prefix) 
            prefix="$2"; 
            shift
            ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 
            ;;
    esac
    shift
done

# Check mutual exclusivity
if [[ $fastq_set -eq 1 && $bam_set -eq 1 ]]; then
    echo "Error: You can only specify either --fq or --bam, not both."
    exit 1
fi

# Check if the correct number of arguments is provided
if { [ -z "$fastq" ] && [ -z "$bam" ]; } || [ -z "$assembly" ] || [ -z "$mapped_folder" ] || [ -z "$thr" ]; then
  echo "Usage: $0 (--fq <fastq_file> | --bam <bam_file>) --fa <assembly> --of <mapped_folder> --t <threads> [--d <value>] [--min_cov <value>]"
  exit 1
fi

# Check if either fastq or bam file exists
if [ ! -z "$fastq" ] && [ ! -f "$fastq" ]; then
  echo "Fastq file does not exist: $fastq"
  exit 1
fi

if [ ! -z "$bam" ] && [ ! -f "$bam" ]; then
  echo "Bam file does not exist: $bam"
  exit 1
fi

# Check if the fasta file exists
if [ ! -f "$assembly" ]; then
  echo "Assembly file does not exist: $assembly"
  exit 1
fi

# Check if the output folder exists and create it if negative
if [ ! -d "$mapped_folder" ]; then
    mkdir -p "$mapped_folder"
fi

# Get the directory path of the main script
current_dir=$(dirname "$(readlink -f "$0")")

# Get filename from the prefix
filename="$prefix"

# Run BWA-MEM to map fastq to fasta
if [ ! -z "$fastq" ]; then
    #filename=$(basename "$fastq" .fastq.gz)
    bwa mem -t "${thr}" "${assembly}" "${fastq}" | samtools view -bS - | samtools sort -o "${mapped_folder}/${filename}.sorted.bam" -
    samtools index "${mapped_folder}/${filename}.sorted.bam"
    echo "${filename} mapped successfully to ${assembly}"
    echo "Extracting low coverage sequences from ${filename}"
    bam="${mapped_folder}/${filename}.sorted.bam"
    bash "$current_dir/scripts/bam2fasta.sh" "${mapped_folder}/${filename}.sorted.bam" "${assembly}" "${min_cov}" "${min_len}" "${d}" "${mapped_folder}" "${prefix}"
else
    #filename=$(basename "${bam%.bam}")
    #filename=$(basename "${filename%.sorted}")
    echo "No fastq file specified. Skipping BWA-MEM mapping."
    samtools index "${bam}"
    echo "Extracting coverage gaps from ${filename}"
    bam="${bam}"
    bash "$current_dir/scripts/bam2fasta.sh" "${bam}" "${assembly}" "${min_cov}" "${min_len}" "${d}" "${mapped_folder}" "${prefix}"
fi

if [[ ! -s "${mapped_folder}/${filename}-GD.fasta" ]]; then
  echo "The file ${mapped_folder}/${filename}-GD.fasta is empty. GenomeDelta was not able to identify any new region in FASTA file. Check if you prepared everything according to the Manual. Try to change the parameters to reduce stringency"
  exit 1
fi

blastn -query "${mapped_folder}/${filename}-GD.fasta" -subject "${mapped_folder}/${filename}-GD.fasta" -out "${mapped_folder}/${filename}-GD.tmp.blast"  -outfmt 6
awk '($12) >= '"${min_bitscore}" "${mapped_folder}/${filename}-GD.tmp.blast" > "${mapped_folder}/${filename}-GD.blast"
rm "${mapped_folder}/${filename}-GD.tmp.blast"
mkdir "${mapped_folder}/${filename}-GD-clusters/"
echo "Finding repetitive gaps..."
python "$current_dir/scripts/blast2clusters.py" "${mapped_folder}/${filename}-GD.blast" "${mapped_folder}/${filename}-GD.fasta" "${mapped_folder}/${filename}-GD-clusters/"
mv "${mapped_folder}/${filename}-GD-clusters/non_rep.fasta" "${mapped_folder}/${filename}-GD-non_rep.fasta"

echo "Extracting consensus sequences of the invaders..."
# Check if there are any .fasta files in the folder
if ! ls "${mapped_folder}/${filename}-GD-clusters/"*.fasta 1> /dev/null 2>&1; then
    echo "No fasta files found in the folder. Zero repetitive sequences found."
    rm -r "${mapped_folder}/${filename}-GD-clusters/"
    exit 1
fi

# Loop through each fasta file in the folder
for fasta in "${mapped_folder}/${filename}-GD-clusters/"*.fasta
do
    # Define the output file names based on the input file name
    output_MSA="${fasta%.fasta}.MSA"
    output_consensus="${fasta%.fasta}.consensus"
    output_standard="${fasta%.fasta}"

    # Run MUSCLE with input and output files
    muscle -in "${output_standard}.fasta" -out "$output_MSA"
    python "$current_dir/scripts/MSA2consensus.py" "$output_MSA" "$output_consensus"
    rm "$output_MSA"
done

# Check if --refine option is activated
if [[ "$refine_set" -eq 1 ]]; then
    echo "Refining..."
    mkdir "${mapped_folder}/${filename}-GD-clusters-refined"
    bash "$current_dir/scripts/find-coupled-clusters.sh" "${mapped_folder}/${filename}-GD-clusters" "${mapped_folder}/${filename}-GD-clusters-refined" "${refine_d}" "${assembly}" "${bam}" "${mapped_folder}/${filename}.bedgraph"
fi

rm "${mapped_folder}/${filename}-flanking.bed"
rm "${mapped_folder}/${filename}.sorted.bam.bai"
rm "${mapped_folder}/${filename}-GD.blast"
rm "${mapped_folder}/${filename}-flanking.bedgraph"
rm "${mapped_folder}/${filename}-GD-flanking.credibility"
rm "${mapped_folder}/${filename}-GD.fasta-e"
mv "${mapped_folder}/${filename}-GD-credibility.bed" "${mapped_folder}/${filename}-GD.bed"
mv "${mapped_folder}/${filename}.fai" "${mapped_folder}/${filename}-GD.fai"
mv "${mapped_folder}/${filename}-GD.blast.sorted" "${mapped_folder}/${filename}-GD.blast"
rm "${mapped_folder}/${filename}.bedgraph"

# Concatenate all consensus files into one candidates file
cat "${mapped_folder}/${filename}-GD-clusters/"*consensus > "${mapped_folder}/${filename}-GD-candidates.fasta"

if [[ "$refine_set" -eq 1 ]]; then
    if [ -n "$(find "${mapped_folder}/${filename}-GD-clusters-refined" -maxdepth 1 -type f -name '*.fasta')" ]; then
        cat "${mapped_folder}/${filename}-GD-clusters-refined/"*consensus >> "${mapped_folder}/${filename}-GD-candidates.fasta"
    fi
fi

# Visualization
samtools faidx "${mapped_folder}/${filename}-GD-candidates.fasta"
samtools faidx "${mapped_folder}/${filename}-GD-non_rep.fasta"
Rscript "$current_dir/scripts/visualization.R" --rep "${mapped_folder}/${filename}-GD-candidates.fasta.fai" --nonrep "${mapped_folder}/${filename}-GD-non_rep.fasta.fai" --output1 "${mapped_folder}/${filename}-GD-candidates.png" --output2 "${mapped_folder}/${filename}-GD-non_rep.png"