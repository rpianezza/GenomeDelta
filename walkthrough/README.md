GenomeDelta - Walkthrough
================

This walkthrough will guide you through an example run of GenomeDelta.

First, install GenomeDelta as specified in the manual.

    cd Applications
    bash ../Downloads/setup.sh

Second, download the input files from the folder
<https://github.com/rpianezza/GenomeDelta/tree/main/walkthrough>:

- test-assembly.fasta
- test-reads.fq

The **fasta** file represents an assembly of a recently collected
specimen, while the **fq** file contains sequencing reads of an old
sample, both from the same species. Now we will run GenomeDelta to
detect if any HT happened between in the species between the two
collection times of our samples.

Activate the conda environment, automatically generated during the
installation.

    conda activate GenomeDelta

Index the fasta file.

    bwa index test-assembly.fasta

Call GenomeDelta using our input files. Please specify:

- –fa: the path to the file “test-assembly.fasta”
- –fq: the path to the file “test-reads.fq”
- –of: the path to folder where you want the GD output
- –prefix: the prefix you want to give to all the GD output files
- –t: the number of threads to use (be careful to use a number of
  threads available in your computer!)

<!-- -->

    GenomeDelta --fa test-assembly.fasta --fq test-reads.fq --of output --prefix test --t 5

You can follow the run looking at the intermediate files being written
in the specified output folder. At the end of the run, your output
folder should contain:

- **.bam** file –\> obtained by mapping the fq file to the fasta file
- **GD.fasta** file –\> the low-coverage sequences extracted from the
  bam file
- **GD.fai** file –\> index file obtained from the GD.fasta
- **GD.blast** file –\> the BLASTn results obtained by aligning all of
  the gaps with each other. Used to identify repetitive clusters.
- **GD.bed** file –\> the low-coverage regions extracted from the bam
  file, in bed format
- **candidates.fasta** and **non_rep.fasta** files –\> repetitive gaps
  and non-repetitive gaps fasta sequences, with their respective
  **.fai** files used for visualization. In this case, the
  **non_rep.fasta** should be empty.
- **.png** file –\> repetitive gaps visualization
- **GD-clusters** folder –\> contains additional information on the
  repetitive clusters found, like the single consensus sequences and the
  fasta file containing every single gap with genomic position and
  coverage bias score.

The visualization (.png) shows a single point (0), with 30 insertions
(y-axis) and 5000 bp (x-axis). The consensus sequence of this cluster is
present in the **candidates.fasta** file (cluster_0), and the single
insertions sequences and genomic locations can be found in the
**GD-clusters** folder.

We found a strong candidate for a recent HTT event.
