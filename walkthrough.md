GenomeDelta - Walkthrough
================

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
theme_set(theme_bw())
```

## How to call GenomeDelta

Example call in the UNIX command line:

    bash main.sh --fq reads.fastq.gz --fa assembly.fa --of folder_path --t 20

The above call is composed of:

- `bash` -\> to specify that we are about to call a bash script
- `main.sh` -\> the main script of **GenomeDelta**.
- `--fq` -\> the **FASTQ** file of the “old” genome.
- `--fa` -\> the **FASTA** file (sometimes only specified as “fa”) of
  the “new” genome (the assembly).
- `--of` -\> the output folder.
- `--t` -\> the number of threads that will be used to parallelize the
  slowest steps (be careful to not use too many threads!).

Remember to index the FASTA assembly before!

    bwa index assembly.fa

### Optional arguments

**GenomeDelta** also has some other options, that can be used to refine
or explore your findings:

- `--min_len` -\> Set the minimum length of a low-coverage region to be
  included in the output files. **Default = 1000**
- `--min_cov` -\> Set the minimum coverage. Below that, a position is
  considered to be “low-coverage” and will be included in the next
  steps. **Default = 2**
- `--d` -\> Set the maximum distance between two low-coverage regions to
  be merged. If the distance between the two regions is below **d**, the
  two regions will be merged. Increasing this distance could create
  artifacts and chimeric sequences, but could find more fragmented
  regions (es. re-invading TEs). **Default = 100**.

### Call GenomeDelta giving a BAM file as input

**GenomeDelta** allows the user to give directly the **BAM** file as
input, instead of the **FASTQ**. These two arguments are mutually
exclusive, so provide only one of the two. Since the mapping is the most
time-consuming step, the user can map once the FASTQ to the FASTA and
subsequently play with the options (–d, –min_len, –min_cov) giving the
already mapped file as input to save time.

Example call:

    bash main.sh --bam reads.fastq.gz --fa assembly.fa --of folder_path --t 5

### Call GenomeDelta giving multiple FASTQ/BAM files as input

To iterate over multiple **FASTQ** or **BAM** files in a folder and run
**GenomeDelta** on all of them against a single assembly, you can use
this one-liner structure:

For **FASTQ**:

    for fq_file in folder/*.fastq.gz; do base_name=$(basename "$fq_file" .fastq.gz); file="$base_name"; bash main.sh --fq "$fq_file" --fa assembly.fa --of folder_path/"$file" --t 20; done

For **BAM** (sorted):

    for bam_file in folder/*.sorted.bam; do base_name=$(basename "$bam_file" .sorted.bam); file="$base_name"; bash main.sh --bam "$bam_file" --fa assembly.fa --of folder_path/"$file" --t 20; done

## Step-by-step explanation of GenomeDelta’s workflow

1)  Mapping the FASTQ old genome to the FASTA recent assembly (tools:
    `bwa-mem`).

2)  Sorting and indexing the mapped file, obtaining the final **bam**
    file (tools: `samtools`).

3)  Calls the script `bam2fasta`, which:

- extracts the low coverage regions of the bam file (coverage \<2,
  meaning only regions with 0 and 1 coverage will be extracted), which
  are the sequences absent from the old genome but present in the new
  one. The output of this step is the **low_coverage.bedgraph** file
  (tools: `samtools`, `awk`).
- merge contiguous bases with low coverage. Then, merge adjacent
  low-coverage regions (gap permitted by default: 100 bases). This is
  done to defragment sequences which may have some reads mapping on a
  fragment in between them. Output: **low_coverage_merged.bed** (tools:
  `bedtools`).
- remove small sequences (\<1000 bases). Output: **unmapped.bed**.
- extract FASTA sequences of the low-coverage regions. Output:
  **unmapped.fasta** (tools: `bedtools`).

4)  BLAST the resulting FASTA file against itself and removes low
    quality hits (BLAST score \< 1000). Output:
    **unmapped-filtered.blast**.

5)  Creates the folder **clusters**.

6)  Calls the script `blast2clusters`, which:

- extracts similar sequences from the **unmapped-filtered.blast** file.
- writes a fasta file for each cluster in the folder **clusters**, but
  only if the cluster contains at least 3 sequences (\>2). This is done
  to keep only repetitive sequences in the final file.

7)  On each **cluster.fasta** file, runs a MSA with `MUSCLE`.

8)  Calls the `MSA2consensus`, which creates a consensus sequence for
    each of the identified clusters.

9)  Concatenate the consensus sequences into the final output with the
    candidates (**candidates.fasta**).

## Output files

- `unmapped.fasta` -\> FASTA containing all the sequences with coverage
  lower than `--min_cov` (default = 2) and longer than `--min_len`
  (default = 1000 bases), merged together if their distance is lower
  than `--d` (default = 100 bases).

- `unmapped.bed` -\> Chromosome, starting and ending positions of the
  sequences collected in `unmapped.fasta`.

- `unmapped.blast` -\> Output of the self BLAST of `unmapped.fasta`
  against itself. This file is then used by the program to find the
  repetitive clusters of sequences.

- `candidates.fasta` -\> Consensus sequences of the repetitive clusters.
  Represents a list of candidates of the invading TEs, to be further
  investigated.

- `clusters` -\> This folder contains, for each of the repetitive
  clusters found:

  - a **FASTA** file containing all the sequences clustering together.
  - a **MSA** (multiple sequence alignment) of the cluster sequences,
    performed using MUSCLE.
  - the **consensus** sequence of the cluster, then concatenated with
    the other consensus into the **candidates.fasta** output.
