GenomeDelta - Manual
================

## Purpose of GenomeDelta

**GenomeDelta** is a software designed to discover horizontal transfer
events in a species. By comparing an older genome (in FASTQ or BAM
format) with a more recent assembly (in FASTA format) of the same
species, **GenomeDelta** identifies novel genetic elements present in
the new assembly, but absent in the older genome.

Its primary focus is on detecting transposable element invasions, while
also offering the capability to unveil other HT events and other genomic
alterations. An enormous advantage of **GenomeDelta** when dealing with
transposable elements is that it does not require any reference library
of transposons to identify the novel invaders.

## Install GenomeDelta

### MacOS

If you are using a MacOS machine, download the `setup.sh` file only
from:
<https://github.com/rpianezza/GenomeDelta/blob/main/macOS/setup.sh>.

Double clicking on the setup file will install GD in the root folder.
Alternatively, open a terminal and move to the folder where you want to
install GD, for example the Applications folder. From the Application
folder, call the downloaded `setup.sh` file:

    cd Applications
    bash ../Downloads/setup.sh

For a **manual installation**, clone the GD repository to your local
folder. Then create the conda environment using the `set-env.yml`.

    git clone https://github.com/rpianezza/GenomeDelta.git
    conda env create -f set-env.yml

### Linux

Download the GD repository using git clone, then create the conda
environment using the `set-env.yml`.

    git clone https://github.com/rpianezza/GenomeDelta.git
    conda-env create -f linux/set-env-linux.yml

## Call GenomeDelta

Activate the conda environment:

    conda activate GenomeDelta

When being into the GenomeDelta conda environment, to run the tool you
can just type GenomeDelta (MacOS) or call the main script defining its
path (if you followed a **manual installation** or you are on a Linux
machine).

Example call (MacOS):

    GenomeDelta --fq reads.fastq.gz --fa assembly.fa --of folder_path --prefix name --t 20

Example call (Linux or manual installation):

    bash main.sh --fq reads.fastq.gz --fa assembly.fa --of folder_path --prefix name --t 20

The above call is composed of:

- `bash` -\> to specify that we are about to call a bash script
- `main.sh` -\> the main script of **GenomeDelta**.
- `--fq` -\> the **FASTQ** file of the “old” genome.
- `--fa` -\> the **FASTA** file of the “new” genome (the assembly).
- `--of` -\> the output folder.
- `--prefix` -\> the prefix that all the output files will have.
- `--t` -\> the number of threads that will be used to parallelize the
  slowest steps.

Remember to index the FASTA assembly before the call!

    bwa index assembly.fa

GD can also accept sorted **BAM** files as input instead of the FASTQ
file. The BAM file should have been mapped to the same FASTA assembly
specified in the call.

    GenomeDelta --bam mapped.sorted.bam --fa assembly.fa --of folder_path --prefix name --t 20

## Optional arguments

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
- `--refine` -\> if the option “refine” is activated, GD will run as
  normal but in the end will try to merge clusters which insertions are
  close to each other. To activate the “refine” option, simply type
  “–refine” at the end of a normal GD call.
- `--refine_d` -\> only if the “refine” option is activated, you can
  specify the maximum distance at which GD will try to match insertions
  in two clusters. **Default = 2500**.

## Output files

- `GD-candidates.png` -\> Visualization of the repetitive clusters found
  (candidate TEs).

- `GD-candidates.fasta` -\> Consensus sequences of the repetitive
  clusters. Represents a list of candidates of the invading TEs. Each
  sequence name is composed off: (i) cluster name, (ii) median
  credibility score of the sequences in the cluster, (iii) number of
  sequences in the cluster.

### Secondary output files

- `GD.fasta` -\> FASTA containing all the sequences with coverage lower
  than `--min_cov` (default = 2) and longer than `--min_len` (default =
  1000 bases), merged together if their distance is lower than `--d`
  (default = 100 bases).

- `GD.bed` -\> Chromosome, starting and ending positions of the
  sequences collected in `GD.fasta`.

- `GD-non_rep.fasta` -\> Subset of the `GD.fasta` file containing only
  sequences which were not included in the repetitive clusters, so
  either non repetitive sequences, huge gaps (\>25.000 bp) or sequences
  which have some similarity with a repetitive clusters but with
  dimensions bigger the the rest of the cluster.

- `GD-non_rep.png` -\> Visualization of the non repetitive gaps found
  (candidate HTs).

- `GD.blast` -\> Output of the self BLAST of `GD.fasta` against itself.
  This file is then used by the program to find the repetitive clusters
  of sequences.

- `GD-clusters` -\> This folder contains, for each of the repetitive
  clusters found:

  - a **FASTA** file containing all the sequences clustering together.
  - the **consensus** sequence of the cluster, then concatenated with
    the other consensus into the **GD-candidates.fasta** file.

## Tips & Tricks

### Call GenomeDelta giving multiple FASTQ/BAM files as input

Example call:

    bash main.sh --bam reads.sorted.bam --fa assembly.fa --of folder_path --prefix name --t 20

### Call GenomeDelta giving multiple FASTQ/BAM files as input

To iterate over multiple **FASTQ** or **BAM** files in a folder and run
**GenomeDelta** on all of them against a single assembly, you can use
this one-liner structure:

For **FASTQ**:

    for fq_file in folder/*.fastq.gz; do base_name=$(basename "$fq_file" .fastq.gz); file="$base_name"; bash main.sh --fq "$fq_file" --fa assembly.fa --of folder_path/"$file" --prefix name --t 20; done

For **BAM** (sorted):

    for bam_file in folder/*.sorted.bam; do base_name=$(basename "$bam_file" .sorted.bam); file="$base_name"; bash main.sh --bam "$bam_file" --fa assembly.fa --of folder_path/"$file" --prefix name --t 20; done

These commands will generate a separate folder for each of the input
files, named as the input file basename.

### Call GenomeDelta giving multiple FASTA assemblies as input

To iterate over multiple **FASTA assemblies** and run **GenomeDelta** on
all of them against a single FASTQ file, you can use this loop
structure:

    for fa_file in /path/to/your/fasta/files/*.fa; do
        base_name=$(basename "$fa_file" .fa)
        file="$base_name"
        bash main.sh --fq reads.fastq.gz --fa "$fa_file" --of /path/to/output/"$file" --prefix name --t 20
    done

This command will generate a separate folder for each of the assemblies,
named as the assembly file basename.

### Call GenomeDelta giving multiple FASTQ/BAM files and multiple FASTA assemblies as input

To iterate over multiple **FASTA assemblies** and run **GenomeDelta** on
all of them against multiple FASTQ files, you can use this double loop
structure:

    for fq_file in /path/to/your/fastq/files/*.fastq.gz; do
        base_name=$(basename "$fq_file" .fastq.gz)
        for fa_file in /path/to/your/fasta/assemblies/*.fa; do
            base_name_fa=$(basename "$fa_file" .fa)
            file="$base_name"
            file_fa="$base_name_fa"
            bash main.sh --fq "$fq_file" --fa "$fa_file" --of /path/to/output/"$file"_"$file_fa" --prefix name --t 20
        done
    done

This command will generate a separate folder for each of the
assembly-FASTQ combination, named as the FASTQ file basename and the
assembly file basename separated by “\_“. You can change this loop to
make it iterate over BAM files instead of FASTQ. Note that you may need
to adjust the extension”fa” to “fasta” based on the assemblies names.

## Interpreting the results and potential issues

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
TE very similar to the old TE, where the reads are mapping). The option
**–refine** is designed to tackle this issue.

Those problems could be solved by playing with the options provided by
the software, but multiple run of **GenomeDelta** could be necessary to
polish the results.

Furthermore, the quality of both the FASTQ file and the FASTA assembly
have an impact on the quality of the results.

- If the **FASTQ** file has a low coverage or/and a high gaps
  percentage, many coverage gaps will be identified by GenomeDelta when
  mapping to the assembly.

- If the **FASTA** assembly has a low quality, for example many contigs,
  there could be contigs where the reads are not mapping at all (es.
  contaminations from other organisms). The quality of the long read
  assembly is thus crucial for a smooth **GenomeDelta** run.

In general, **clusters with high credibility scores (close to 1), with a
long consensus sequence and composed by many sequences are likely to be
more valuable.**

## Covarage bias score

Each extracted sequence has its own coverage bias score, calculated as
the proportion between the coverage of the 10kb flanking regions and the
overall mean coverage. The minimum score is 0, while a credibility score
close to 1 or above 1 is highly rated. The coverage bias scores are
included in the FASTA name of the sequences, which is in the format
**chr:start-end-credibility**. The coverage bias score assigned to the
repetitive clusters is the median coverage bias scores of the sequences
in the cluster.
