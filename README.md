# GenomicUF
Scripts to implement Genomic UniFrac for Celeste's adaptation project. Eventually will scale so that code is more general.

## What is Genomic UniFrac?
Genomic UniFrac implements the [UniFrac algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1317376/) on whole genome sequences. More than one microbial species may be sequenced (i.e. the data may be metagenomic in neature), but whole genome sequences are available. This is distinct from Meta UniFrac (or Metagenomic UniFrac) which operates on sequencing data obtained from a mixed microbial community but whole genome sequences are not necessarily inferred. 

## Required inputs
* .fasta / .fa / .fna files containing genome sequences 
* .gff files for each genome sequence

## General analysis pipeline
* Quality control on GFF (general feature format; tab delimited file containing genomic features) files obtained from sequences
* Concatenate GFF files into a pandas DataFrame for ease of parsing information
* Align genomic sequences using [HMMER3](http://hmmer.org/).
* Perform multiple sequence alignment
* 
