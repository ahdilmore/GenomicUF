# GenomicUF
Scripts to process data for [meta]genomic UniFrac calculations. 

## What is Genomic UniFrac?
Genomic UniFrac implements the [UniFrac algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1317376/) on whole genome sequences. More than one microbial species may be sequenced (i.e. the data may be metagenomic in neature), but whole genome sequences are available. This is distinct from Meta UniFrac (or Metagenomic UniFrac) which operates on sequencing data obtained from a mixed microbial community but whole genome sequences are not necessarily inferred. 

## Required inputs (for full processing)
* .fasta / .fa / .fna files containing genome sequences 
* .gff (general feature format; tab delimited file containing genomic features) files for each genome sequence
