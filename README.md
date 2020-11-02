# Genomic_population_analysis_of_Neisseria_gonorrhoeae
This repository contains the Snakemake pipeline used to generate the results described in the paper: "Emergence of a cephalosporin reduced susceptible Neisseria gonorrhoeae clone between 2014-2019 in Amsterdam, the Netherlands, revealed by a genomic population analysis", preprint:

## Input
This pipeline uses forward and reverse raw Illumina sequencing reads which are located in the folder `raw_data`. The raw data files should be named `{id}_R1.fastq.gz` and `{id}_R2.fastq.gz`

## Pipeline 
The pipeline includes the following steps and tools:

| Step     | Tool     |   
| ---------|----------|
| Filter low quality raw reads + trim adapters | fastp |
| Assembly | skesa | 
| Assembly quality check | QUAST |
| Calculate percentage of reference FA1090 covered by raw reads | minimap2 + samtools |
| Multi-Locus Sequence Typing | mlst |
| NG-Multi Antigen Sequence Typing | ngmaster |
| Look for23S mutations in raw reads | ariba |
| Call variants using reference FA1090 + create core genome alignment | snippy |
| Remove recombination from variant alignment and create phylogenetic tree | Gubbins + RAxML |

## Other files needed
- Reference genome FA1090 is used for calculating coverage and calling variants: NC_002946.2. The reference genome should be located in the same directory as the Snakefile.
- Reference 23S database: `23S_seq.fa` used to create database with `ariba` using command `ariba prepareref --all_coding no -f 23S_seq.fa ariba/23S_seq.prepareref`. For running the pipeline the reference database should be located as `ariba/23S_seq.prepareref`.

