GenomeDelta - Manual
================

# Purpose of GenomeDelta

**GenomeDelta** is a software designed to unravel the mysteries of
genome evolution in a species. By comparing an older genome (in FASTQ
format) with a more recent assembly (in FASTA format) of the same
species, **GenomeDelta** swiftly identifies novel genetic elements
present in the new assembly, but absent in the older genome.

Its primary focus is on detecting transposable element invasions, while
also offering the capability to pinpoint various other noteworthy
genomic alterations. An enormous advantage of **GenomeDelta** when
dealing with transposable elements is that it does not require any
reference library of transposons to identify the novel invaders.
Reference libraries are often incomplete and hard to build, especially
for non-model species.

# Set the GenomeDelta conda environment

First, you need to create the **conda environment** to install all the
necessary packages. Use the `set-env.yml` file.

If you are using a MacOS machine:

    conda env create -f macOS/set-env-MAC.yml

If you are using a Linux machine:

    conda-env create -f linux/set-env-linux.yml

Then, before calling GenomeDelta, activate the environment:

    conda activate GenomeDelta

# Call GenomeDelta

After you are sure to be into the GenomeDelta conda environment, you can
call the main script. The main script is located in the folder **Linux**
or **macOS**, use the one you need. Example call in the UNIX command
line:

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
  steps. **Default = 1**
- `--d` -\> Set the maximum distance between two low-coverage regions to
  be merged. If the distance between the two regions is below **d**, the
  two regions will be merged. Increasing this distance could create
  artifacts and chimeric sequences, but could find more fragmented
  regions (es. re-invading TEs). **Default = 100**.
- `--min_bitscore` -\> to find repetitive clusters, GD is using BLAST.
  The output is filtered based on the **bitscore** value, with default
  set to 1000 to only consider high quality alignments. If you want to
  find small sequences, you may want to decrease this parameter.

### Call GenomeDelta giving a BAM file as input

**GenomeDelta** allows the user to give directly the **BAM** file as
input, instead of the **FASTQ**. These two arguments are mutually
exclusive, so provide only one of the two. Since the mapping is the most
time-consuming step, the user can map once the FASTQ to the FASTA and
subsequently play with the options (–d, –min_len, –min_cov) giving the
already mapped file as input to save time.

Example call:

    bash main.sh --bam reads.sorted.bam --fa assembly.fa --of folder_path --t 20

### Call GenomeDelta giving multiple FASTQ/BAM files as input

To iterate over multiple **FASTQ** or **BAM** files in a folder and run
**GenomeDelta** on all of them against a single assembly, you can use
this one-liner structure:

For **FASTQ**:

    for fq_file in folder/*.fastq.gz; do base_name=$(basename "$fq_file" .fastq.gz); file="$base_name"; bash main.sh --fq "$fq_file" --fa assembly.fa --of folder_path/"$file" --t 20; done

For **BAM** (sorted):

    for bam_file in folder/*.sorted.bam; do base_name=$(basename "$bam_file" .sorted.bam); file="$base_name"; bash main.sh --bam "$bam_file" --fa assembly.fa --of folder_path/"$file" --t 20; done

These commands will generate a separate folder for each of the input
files, named as the input file basename.

### Call GenomeDelta giving multiple FASTA assemblies as input

To iterate over multiple **FASTA assemblies** and run **GenomeDelta** on
all of them against a single FASTQ file, you can use this loop
structure:

    for fa_file in /path/to/your/fasta/files/*.fa; do
        base_name=$(basename "$fa_file" .fa)
        file="$base_name"
        bash main.sh --fq reads.fastq.gz --fa "$fa_file" --of /path/to/output/"$file" --t 20
    done

This command will generate a separate folder for each of the assemblies,
named as the assembly file basename.

### Call GenomeDelta giving multiple FASTQ/BAM files as well as multiple FASTA assemblies as input

To iterate over multiple **FASTA assemblies** and run **GenomeDelta** on
all of them against multiple FASTQ files, you can use this double loop
structure:

    for fq_file in /path/to/your/fastq/files/*.fastq.gz; do
        base_name=$(basename "$fq_file" .fastq.gz)
        for fa_file in /path/to/your/fasta/assemblies/*.fa; do
            base_name_fa=$(basename "$fa_file" .fa)
            file="$base_name"
            file_fa="$base_name_fa"
            bash main.sh --fq "$fq_file" --fa "$fa_file" --of /path/to/output/"$file"_"$file_fa" --t 20
        done
    done

This command will generate a separate folder for each of the
assembly-FASTQ combination, named as the FASTQ file basename and the
assembly file basename separated by “\_“. You can change this loop to
make it iterate over BAM files instead of FASTQ. Note that you may need
to adjust the extension”fa” to “fasta” based on the assemblies names.

# Output files

- `GD.fasta` -\> FASTA containing all the sequences with coverage lower
  than `--min_cov` (default = 2) and longer than `--min_len` (default =
  1000 bases), merged together if their distance is lower than `--d`
  (default = 100 bases).

- `GD-non_rep.fasta` -\> Subset of the `GD.fasta` file containing only
  sequences which were not included in the repetitive clusters, so
  either non repetitive sequences, huge gaps (\>25.000 bp) or sequences
  which have some similarity with a repetitive clusters but with
  dimensions bigger the the rest of the cluster.

- `GD.bed` -\> Chromosome, starting and ending positions of the
  sequences collected in `GD.fasta`.

- `GD.blast` -\> Output of the self BLAST of `GD.fasta` against itself.
  This file is then used by the program to find the repetitive clusters
  of sequences.

- `GD-candidates.fasta` -\> Consensus sequences of the repetitive
  clusters. Represents a list of candidates of the invading TEs, to be
  further investigated.

- `GD-clusters` -\> This folder contains, for each of the repetitive
  clusters found:

  - a **FASTA** file containing all the sequences clustering together.
  - a **MSA** (multiple sequence alignment) of the cluster sequences,
    performed using MUSCLE.
  - the **consensus** sequence of the cluster, then concatenated with
    the other consensus into the **candidates.fasta** output.

# Interpreting the results and potential issues

The main output file for repetitive sequences is `GD-candidates.fasta`,
which represents a list of candidates invaders in the time range going
from the old genome collection time to the recent assembly collection
time. However, not all the sequences found in the file are probably real
invaders.

Inside the file, you will probably find sequences of **telomeric
repeats**: telomeres often have coverage gaps in sequencing efforts,
thus a telomeric repetitive region could remain unsequenced in the FASTQ
file and appear as coverage gap in the FASTA assembly, ending in our
candidates list. It is thus important to check where are the repetitive
sequences of the cluster (the FASTA name of the sequences contains the
genomic position).

Also, an invading TE with a consistent sequence similarity to another
(old, non recently invading) TE will result in reads mapping to the
supposed “coverage gap” of the new TE. Thus, instead of having a clean
coverage gap representing the new TE insertion, we will have some reads
mapping on at least a part of the gap. In the **GenomeDelta** output,
you may find **different parts of the same, new TE in different
“clusters”**, so in different FASTA entries, each representing a “clean”
gap in the coverage divided by a high coverage region (es. a part of the
TE very similar to the old TE, where the reads are mapping).

Those problems could be solved by playing with the options provided by
the software, but multiple run of **GenomeDelta** could be necessary to
polish the results.

Furthermore, the quality of both the FASTQ file and the FASTA assembly
have an impact on the purity of the results.

- If the **FASTQ** file has a low coverage or/and a high gaps
  percentage, many coverage gaps will be identified by GenomeDelta when
  mapping to the assembly.

- If the **FASTA** assembly has a low quality, for example many contigs,
  there could be contigs where the reads are not mapping at all (es.
  contaminations from other organisms). The quality of the long read
  assembly is thus crucial for a smooth **GenomeDelta** run.

## Credibility scores

Each extracted sequence has its own credibility score, calculated as the
proportion between the coverage of the 10kb flanking regions and the
overall mean coverage. The minimum score is 0, while a credibility score
close to 1 or above 1 is highly rated. The credibility scores are
included in the FASTA name of the sequences, which is in the format
**chr:start-end-credibility**. The credibility score assigned to the
repetitive clusters is the average credibility scores of the sequences
in the cluster.

# Step-by-step explanation of GenomeDelta’s workflow

1)  Mapping the FASTQ old genome to the FASTA recent assembly (tools:
    `bwa-mem`).

2)  Sorting and indexing the mapped file, obtaining the final **bam**
    file (tools: `samtools`).

3)  Calls the script `bam2fasta`, which:

- extracts the low coverage regions of the bam file (with default
  coverage \<1, meaning only regions with 0 coverage will be extracted),
  which are the sequences absent from the old genome but present in the
  new one. The output of this step is the **low_coverage.bedgraph** file
  (tools: `samtools`, `awk`).
- merge contiguous bases with low coverage. Then, merge adjacent
  low-coverage regions (gap permitted by default: 100 bases). This is
  done to defragment sequences which may have some reads mapping on a
  fragment in between them. Output: **low_coverage_merged.bed** (tools:
  `bedtools`).
- remove small sequences (\<1000 bases). Output: **GD.bed**.
- extract FASTA sequences of the low-coverage regions. Output:
  **GD.fasta** (tools: `bedtools`).

4)  BLAST the resulting FASTA file against itself and removes low
    quality hits (BLAST score \< 1000). Output: **GD.blast**.

5)  Creates the folder **GD-clusters**.

6)  Calls the script `blast2clusters`, which:

- extracts similar sequences from the **GD.blast** file.
- writes a fasta file for each cluster in the folder **GD-clusters**,
  but only if the cluster contains at least 3 sequences (\>2). This is
  done to keep only repetitive sequences in the final file.

7)  On each **cluster.fasta** file, runs a MSA with `MUSCLE`.

8)  Calls the `MSA2consensus`, which creates a consensus sequence for
    each of the identified clusters using a majority-wins approach.

9)  Concatenate the consensus sequences into the final output with the
    candidates (**GD-candidates.fasta**).

10) All the extracted sequences which are not included in the clusters
    can be found in the **GD-non_rep.fasta**
