Testing Scenarios
================
2023-11-03

## Overview

GenomeDelta is able to identify transposable element invasions by
mapping the genome assembly of interest (fasta format) to a reference
genome of the same species (fastq format) and extracting the gaps found
in the assembly. For the purpose of testing GenomeDelta, several
artificial scenarios of Transposon invasions were created using simulaTE
v1.13.

### Artificial Data

The testing was performed with Drosophila melanogaster Chromosome 1 as
reference assembly and Transposon sequences of size 5000 bp, 2500 bp,
1000 bp, 500 bp and 250 bp, that were randomly generated using the
Random DNA Sequence Generator service provided by the Riverside
University of California
(<http://www.faculty.ucr.edu/~mmaduro/random.htm>).

### Read generation-Standard Parameters

Templatereads were generated using `create-reads.py` with the following
standard parameters:

``` shell
python2 create-reads.py --fasta assembly.fa --coverage 10 --read-length 100 --output templatereads.fq --method random [--error-rate 0]``
```

Resulting in an output file containing 100 bp long randomly distributed
reads of the chasis with a mean coverage of 10 reads per 100 bp and
errror-rate 0.

Below the definition of each parameter:

- `--fasta` requires a template genome in .fasta/.fa format

- `--coverage` defining the mean coverage throughout the template genome

- `--read-length` The length of each read in bp

- `--output` Name of the ouput file

- `--method` Choose if the created reads are randomly distributed or
  uniform. Can be either `random`or `uniform`.

  Optional argument

- `--error-rate` A value that can determine the frequency of erroneous
  base-implementaion.

### Using simulaTE to create artificial invasions

To generate the artificial TE invasions necessary to test GenomeDelta,
two scripts contained in simulaTE 1.13, namely
`define-landscape_random-insertions-freq-range.py` and
`build-population-genome.py` were used as follows:

``` shell
python2 define-landscape_random-insertions-freq-range.py --assembly.fa --te-seqs dummyTE.fasta --insert-count 25 --min-freq 1 
--max-freq 1 --min-distance 500 --N 1 --output dummyTE_invasion.pgd
```

and

``` shell
python2 build-population-genome.py --pgd dummyTE_invasion.pgd --te-seqs dummyTE.fasta --chassis assembly.fa --output dummyTE_invasion.fasta
```

### Using GenomeDelta to detect the artificial TE insertions

Before running GenomeDelta on an assembly, the indexing of the fasta
file is required.

``` shell
bwa index assembly.fa
```

To run the software, the following command was used:

``` shell
GenomeDelta --fq reads.fq --fa assembly.fa --of folder_path --t 20 --prefix Test [--min_len MINIMUM LENGTH OF GAP --d MAXIMUM DISTANCE BETWEEN GAPS --min_cov MINIMUM READ COVERAGE] 
```

Below the definition of each parameter:

- `--fq` requires the reference genome reads in .fq format

- `--fa` requires the genome of interest assembly in .fa/fasta format

- `--of` defining the output folder `--t` number of threads to be used

- `--prefix` defining the prefix of the output file

- `--t`number of threads that will be used during the mapping step

Optional arguments

- `--min_length` Set the minimum length of a low-coverage region to be
  included in the output files. Default=1000

- `--min_cov` Set the minimum coverage. Below that, a position is
  considered to be “low-coverage” and will be included in the next
  steps. Default=1

- `--d` Set the maximum distance between two low-coverage regions to be
  merged. If the distance between the two regions is below d, the two
  regions will be merged. Increasing this distance could create
  artifacts and chimeric sequences, but could find more fragmented
  regions (es. re-invading TEs). Default = 100.

### Scenario 1 - Random vs Uniform vs aDNA coverage with means 1/5/10

Using both the `random` and `uniform`method, reads with `1`, `5` and
`10` coverage where created, resulting in a total of 6 .fq files. Those
files will take on the role of the “old Genome”. For the new Genome, a
randomly generated dummy transposon of `5000 bp` length was inserted
into the chasis with simulaTE using the above defined standard
parameters.

To simulate ancient DNA reads contaminated with modern Drososphila DNA
as well as DNA from the endosymbiont Wolbachia, we used Gargammel v3.
Below an example to create reads with mean coverage 5, where `-c`
corresponds to the desired coverage, `-comp` gives the DNA composition
of bacterial contaminant, modern contaminant and DNA of interest
respectively and `-o` specifies the path to the ouputfolder.

``` shell
gargammel.pl -c 5 --comp 0.1,08,0.82 -o outputfolder inputfolder
```

### Scenario 2 - Error-rate testing

To check if Errors in the reference genome have a significant effect on
the performance of GenomeDelta, reads with error-rate `0.01` and `0.05`
were created for `coverages 1/5/10` with the `random` method, resulting
in a total of 6 .fq files. Below the example for creating reads with
0.01 error-rate and 1 coverage.

``` shell
python2 create-reads.py --fasta assembly.fa --coverage 1 --read-length 100 --output templatereads.fq --method random --error-rate 0.01
```

### Scenario 3 - Read-length testing

To check if the read-length in the reference genome has a significant
effect on the performance of GenomeDelta, reads with a length of `50`
and `200` were created for `coverages 1/5/10` with the `random` method,
resulting in a total of 6 .fq files. Below the example for creating
reads with 50 read-length and 1 coverage.

``` shell
python2 create-reads.py --fasta assembly.fa --coverage 1 --read-length 50 --output templatereads.fq --method random 
```

### Scenario 4 - Insertion number testing

To check if the number of TE insertions in the genome of interest has a
significant effect on the performance of GenomeDelta, Genomes with `5`
and `10` insertions were created additionally to the standard scenario
of `25` insertions, resulting in three different genomes of interest.
All of them were mapped to `1/5/10 coverage` reference genomes for a
total of nine iterations of scenario 4. Below the example for creating a
genome of interest with 5 insertions:

``` shell
python2 define-landscape_random-insertions-freq-range.py --chassis assembly.fa --te-seqs dummyTE.fasta --insert-count 5 --min-freq 1 
--max-freq 1 --min-distance 500 --N 1 --output dummyTE_invasion_x5.pgd

python2 build-population-genome.py --pgd dummyTE_invasion_x5.pgd --te-seqs dummyTE.fasta --chassis assembly.fa --output dummyTE_invasion_x5.fasta
```

### Scenario 5 - Insertion length testing

To check if the length of TE insertions in the genome of interest has a
significant effect on the performance of GenomeDelta, Genomes with TE
length of `250 bp`, `500 bp`, `1000 bp` and `2500 bp` were created
additionally to the standard scenario of `5000 kb` insertion length,
resulting in five different genomes of interest. All of them were mapped
to `1/5/10 coverage` reference genomes for a total of fifteen iterations
of scenario 5.

### Scenario 6 - Testing the maximum distance between gaps

To check if changing the maximum distance between gaps (Parameter = d)
in the alignment of the assembly to the reference genome affects the
performance of GenomeDelta, iterations with d = `15/50/100` were
compared for insertion lengths in the assembly of `250 bp`, `500 bp`,
`1000 bp`, `2500` `bp` and `5000 bp` and read coverage of 1 in the
reference genome. Below the example for d = 15

``` shell
GenomeDelta --fq reads.fq --fa assembly.fasta --of folder_path --t 20 --d 15
```
