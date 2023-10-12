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

    bash main.sh /Volumes/Storage/data/Dmel-time-series-data/SRR23876562.fastq.gz /Volumes/EXT-RICCARDO/GenomeDelta/test/ref_Es_Ten/D.mel.Es_Ten.fa /Volumes/EXT-RICCARDO/GenomeDelta/test2 5

It is crucial to be located in the folder where the program is located,
with the `main.sh` script and the folder `scripts` containing the
scripts called from the main (use `cd path/to/folder` to go to the right
folder before calling the program).

The above call is composed of:

- `bash` -\> to specify that we are about to call a bash script
- `main.sh` -\> the name of the script to call
- `.../SRR23876562.fastq.gz` -\> the **FASTQ** file of the “old” genome
- `.../D.mel.Es_Ten.fa` -\> the **FASTA** file (sometimes only specified
  as “fa”) of the “new” genome (the assembly)
- `.../test2` -\> the output folder
- `5` -\> the number of threads that will be used to parallelize the
  slowest steps (be careful to not use too many threads!)

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
